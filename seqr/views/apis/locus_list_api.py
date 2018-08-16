import json
import logging

from django.contrib.auth.decorators import login_required
from django.db.models import Q
from django.views.decorators.csrf import csrf_exempt
from guardian.shortcuts import assign_perm, remove_perm, get_objects_for_group


from reference_data.models import GENOME_VERSION_GRCh37
from seqr.models import LocusList, LocusListGene, LocusListInterval, IS_OWNER, CAN_VIEW, CAN_EDIT
from seqr.model_utils import create_seqr_model, delete_seqr_model, find_matching_xbrowse_model
from seqr.utils.xpos_utils import get_xpos
from seqr.views.apis.auth_api import API_LOGIN_REQUIRED_URL
from seqr.views.utils.gene_utils import get_genes, get_gene_symbols_to_gene_ids
from seqr.views.utils.json_utils import create_json_response
from seqr.views.utils.json_to_orm_utils import update_model_from_json
from seqr.views.utils.orm_to_json_utils import get_json_for_locus_lists, get_json_for_locus_list
from seqr.views.utils.permissions_utils import get_project_and_check_permissions, check_object_permissions, check_public_object_permissions
from xbrowse_server.base.models import ProjectGeneList as BaseProjectGeneList

logger = logging.getLogger(__name__)

INVALID_ITEMS_ERROR = 'This list contains invalid genes/ intervals. Update them, or select the "Ignore invalid genes and intervals" checkbox to ignore.'


@login_required(login_url=API_LOGIN_REQUIRED_URL)
@csrf_exempt
def locus_lists(request):
    locus_lists = LocusList.objects.filter(Q(is_public=True) | Q(created_by=request.user))
    locus_lists_json = get_json_for_locus_lists(locus_lists, request.user)

    return create_json_response({
        'locusListsByGuid': {locus_list['locusListGuid']: locus_list for locus_list in locus_lists_json}
    })


@login_required(login_url=API_LOGIN_REQUIRED_URL)
@csrf_exempt
def locus_list_info(request, locus_list_guid):
    locus_list = LocusList.objects.get(guid=locus_list_guid)

    check_public_object_permissions(locus_list, request.user)

    locus_list_json = get_json_for_locus_list(locus_list, request.user)
    gene_ids = [item['geneId'] for item in locus_list_json['items'] if item.get('geneId')]
    return create_json_response({
        'locusListsByGuid': {locus_list_guid: locus_list_json},
        'genesById': get_genes(gene_ids)
    })


@login_required(login_url=API_LOGIN_REQUIRED_URL)
@csrf_exempt
def create_locus_list_handler(request):
    request_json = json.loads(request.body)

    if not request_json.get('name'):
        return create_json_response({}, status=400, reason='"Name" is required')

    new_genes, existing_gene_ids, new_intervals, existing_interval_guids, invalid_items = _parse_list_items(request_json)
    if invalid_items and not request_json.get('ignoreInvalidItems'):
        return create_json_response({'invalidLocusListItems': invalid_items}, status=400, reason=INVALID_ITEMS_ERROR)

    locus_list = create_seqr_model(
        LocusList,
        name=request_json['name'],
        description=request_json.get('description') or '',
        is_public=request_json.get('isPublic') or False,
        created_by=request.user,
    )
    _update_locus_list_items(locus_list, new_genes, existing_gene_ids, new_intervals, existing_interval_guids,
                             request_json, request.user)
    add_locus_list_user_permissions(locus_list)

    return create_json_response({
        'locusListsByGuid': {locus_list.guid: get_json_for_locus_list(locus_list, request.user)},
        'genesById': new_genes,
    })


@login_required(login_url=API_LOGIN_REQUIRED_URL)
@csrf_exempt
def update_locus_list_handler(request, locus_list_guid):
    locus_list = LocusList.objects.get(guid=locus_list_guid)
    check_object_permissions(locus_list, request.user, permission_level=CAN_EDIT)

    request_json = json.loads(request.body)

    new_genes, existing_gene_ids, new_intervals, existing_interval_guids, invalid_items = _parse_list_items(request_json)
    if invalid_items and not request_json.get('ignoreInvalidItems'):
        return create_json_response({'invalidLocusListItems': invalid_items}, status=400, reason=INVALID_ITEMS_ERROR)

    update_model_from_json(locus_list, request_json, allow_unknown_keys=True)
    _update_locus_list_items(locus_list, new_genes, existing_gene_ids, new_intervals, existing_interval_guids,
                             request_json, request.user)

    return create_json_response({
        'locusListsByGuid': {locus_list.guid: get_json_for_locus_list(locus_list, request.user)},
        'genesById': new_genes,
    })


@login_required(login_url=API_LOGIN_REQUIRED_URL)
@csrf_exempt
def delete_locus_list_handler(request, locus_list_guid):
    locus_list = LocusList.objects.get(guid=locus_list_guid)
    check_object_permissions(locus_list, request.user, permission_level=CAN_EDIT)

    delete_seqr_model(locus_list)
    return create_json_response({'locusListsByGuid': {locus_list_guid: None}})


@login_required(login_url=API_LOGIN_REQUIRED_URL)
@csrf_exempt
def add_project_locus_lists(request, project_guid):
    project = get_project_and_check_permissions(project_guid, request.user, CAN_EDIT)
    xbrowse_project = find_matching_xbrowse_model(project)
    request_json = json.loads(request.body)
    locus_lists = LocusList.objects.filter(guid__in=request_json['locusListGuids'])
    for locus_list in locus_lists:
        assign_perm(user_or_group=project.can_view_group, perm=CAN_VIEW, obj=locus_list)
        xbrowse_gene_list = find_matching_xbrowse_model(locus_list)
        if xbrowse_project and xbrowse_gene_list:
            BaseProjectGeneList.objects.get_or_create(project=xbrowse_project, gene_list=xbrowse_gene_list)

    return create_json_response({
        'locusListGuids': [locus_list['locusListGuid'] for locus_list in get_sorted_project_locus_lists(project, request.user)],
    })


@login_required(login_url=API_LOGIN_REQUIRED_URL)
@csrf_exempt
def delete_project_locus_lists(request, project_guid):
    project = get_project_and_check_permissions(project_guid, request.user, CAN_EDIT)
    xbrowse_project = find_matching_xbrowse_model(project)
    request_json = json.loads(request.body)
    locus_lists = LocusList.objects.filter(guid__in=request_json['locusListGuids'])
    for locus_list in locus_lists:
        remove_perm(user_or_group=project.can_view_group, perm=CAN_VIEW, obj=locus_list)
        xbrowse_gene_list = find_matching_xbrowse_model(locus_list)
        if xbrowse_project and xbrowse_gene_list:
            BaseProjectGeneList.objects.filter(project=xbrowse_project, gene_list=xbrowse_gene_list).delete()

    return create_json_response({
        'locusListGuids': [locus_list['locusListGuid'] for locus_list in get_sorted_project_locus_lists(project, request.user)],
    })


def get_project_locus_list_models(project):
    return get_objects_for_group(project.can_view_group, CAN_VIEW, LocusList)


def get_sorted_project_locus_lists(project, user):
    result = get_json_for_locus_lists(get_project_locus_list_models(project), user)
    return sorted(result, key=lambda locus_list: locus_list['name'])


def add_locus_list_user_permissions(locus_list):
    assign_perm(user_or_group=locus_list.created_by, perm=IS_OWNER, obj=locus_list)
    assign_perm(user_or_group=locus_list.created_by, perm=CAN_EDIT, obj=locus_list)
    assign_perm(user_or_group=locus_list.created_by, perm=CAN_VIEW, obj=locus_list)


def _parse_list_items(request_json):
    requested_items = (request_json.get('parsedItems') or {}).get('items') or []

    existing_gene_ids = set()
    new_gene_symbols = set()
    new_gene_ids = set()
    existing_interval_guids = set()
    new_intervals = []
    invalid_items = []
    for item in requested_items:
        if item.get('locusListIntervalGuid'):
            existing_interval_guids.add(item.get('locusListIntervalGuid'))
        elif item.get('geneId'):
            if item.get('symbol'):
                existing_gene_ids.add(item.get('geneId'))
            else:
                new_gene_ids.add(item.get('geneId'))
        elif item.get('symbol'):
            new_gene_symbols.add(item.get('symbol'))
        else:
            try:
                item['start'] = int(item['start'])
                item['end'] = int(item['end'])
                if item['start'] > item['end']:
                    raise ValueError
                get_xpos(item['chrom'], int(item['start']))
                new_intervals.append(item)
            except (KeyError, ValueError):
                invalid_items.append('chr{chrom}:{start}-{end}'.format(
                    chrom=item.get('chrom', '?'), start=item.get('start', '?'), end=item.get('end', '?')
                ))

    gene_symbols_to_ids = get_gene_symbols_to_gene_ids(new_gene_symbols)
    invalid_items += [symbol for symbol, gene_id in gene_symbols_to_ids.items() if not gene_id]
    new_genes = get_genes([gene_id for gene_id in gene_symbols_to_ids.values() if gene_id] + list(new_gene_ids))
    invalid_items += [gene_id for gene_id, gene in new_genes.items() if not gene]
    new_genes = {gene_id: gene for gene_id, gene in new_genes.items() if gene}
    return new_genes, existing_gene_ids, new_intervals, existing_interval_guids, invalid_items


def _update_locus_list_items(locus_list, new_genes, existing_gene_ids, new_intervals, existing_interval_guids, request_json, user):
    # Update genes
    for locus_list_gene in locus_list.locuslistgene_set.exclude(gene_id__in=existing_gene_ids):
        delete_seqr_model(locus_list_gene)

    for gene_id in new_genes.keys():
        create_seqr_model(
            LocusListGene,
            locus_list=locus_list,
            gene_id=gene_id,
            created_by=user,
        )

    # Update intervals
    genome_version = request_json.get('intervalGenomeVersion') or GENOME_VERSION_GRCh37
    locus_list.locuslistinterval_set.exclude(guid__in=existing_interval_guids).delete()
    for existing_interval in locus_list.locuslistinterval_set.all():
        update_model_from_json(existing_interval, {'genomeVersion': genome_version})
    for interval in new_intervals:
        LocusListInterval.objects.create(
            locus_list=locus_list,
            chrom=interval['chrom'].lstrip('chr'),
            start=interval['start'],
            end=interval['end'],
            genome_version=genome_version,
        )