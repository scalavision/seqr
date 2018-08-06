import logging

from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.core.management.base import BaseCommand
from django.db.models.query_utils import Q
from tqdm import tqdm

from seqr.models import Project, Family, Sample
from seqr.views.utils.permissions_utils import check_permissions
from xbrowse.core.constants import CLINVAR_ANNOTATION_DEFINITIONS

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = 'Variant search'

    def add_arguments(self, parser):

        parser.add_argument("--user", help="user email. Run the search as this user (for testing permissions)")

        # select families
        parser.add_argument("projects", nargs="*", help="Project name or guid. If not specified, defaults to all projects.")
        parser.add_argument("--project-category", nargs="*", help="Project categories")
        parser.add_argument("--analysis-group", nargs="*", help="Analysis group name or guid")
        parser.add_argument("-f", "--family", nargs="*", help="Family guid")
        parser.add_argument("-i", "--individual", nargs="*", help="Individual guid")

        # inheritance filters
        parser.add_argument("--inheritance-mode", nargs="*", choices=["HOM-REF", "XLINKED", "DOM", "DENOVO"])
        parser.add_argument("--allow-missing-genotypes", action="store_true")

        # QC filters
        parser.add_argument("--pass", action="store_true")
        parser.add_argument("--AB-affected", type=float)
        parser.add_argument("--DP-affected", type=int)
        parser.add_argument("--GQ-affected", type=int)
        # parser.add_argument("--mq", type=float)
        # parser.add_argument("--qd", type=float)

        # loci
        parser.add_argument("--gene-list", nargs="*", help="gene list")
        parser.add_argument("--gene-id", nargs="*", help="gene id")
        parser.add_argument("--transcript-id", nargs="*", help="transcript id")
        parser.add_argument("--interval", nargs="*", help="interval")
        parser.add_argument("--variant", nargs="*", help="chrom-pos-ref-alt")

        # dataset type
        parser.add_argument("--dataset-type", nargs="*", choices=["VARIANT", "SV"])
        parser.add_argument("--sample-type", nargs="*", choices=["WES", "WGS"])
        parser.add_argument("--variant-type", nargs="*", choices=["SNP", "INS", "DEL"])

        parser.add_argument("--dataset-type2", nargs="*", choices=["VARIANT", "SV"])
        parser.add_argument("--sample-type2", nargs="*", choices=["WES", "WGS"])
        parser.add_argument("--variant-type2", nargs="*", choices=["SNP", "INDEL"])

        # deleteriousness predictors
        #  dbnsfp_DANN_score               0.57
        #  dbnsfp_FATHMM_pred               D;D;D;D;D;D;D;D
        #  dbnsfp_GERP_RS               1.87
        #  dbnsfp_MetaSVM_pred               T
        #  dbnsfp_MutationTaster_pred               P;P;P;P;P;P
        #  dbnsfp_Polyphen2_HVAR_pred               B;B;B;B;B
        #  dbnsfp_REVEL_score               0.368
        #  dbnsfp_SIFT_pred               T;T;T;T;T;T;T;T;T
        #  dbnsfp_phastCons100way_vertebrate               0.002000
        #  mpc_MPC               0.811105209025
        #  cadd_PHRED               0.355
        #  eigen_Eigen_phred               0.547
        #  primate_ai_score            ...

        # allele frequencies
        parser.add_argument("--this-callset-AC", type=int)
        parser.add_argument("--this-callset-AF", type=float)

        parser.add_argument("--g1k-AC", type=int)
        parser.add_argument("--g1k-AF", type=float)

        parser.add_argument("--exac-AC", type=int)
        parser.add_argument("--exac-HH", type=int)
        parser.add_argument("--exac-AF", type=float)

        parser.add_argument("--gnomad-exomes-AC", type=int)
        parser.add_argument("--gnomad-exomes-HH", type=int)
        parser.add_argument("--gnomad-exomes-AF", type=float)

        parser.add_argument("--gnomad-genomes-AC", type=int)
        parser.add_argument("--gnomad-genomes-HH", type=int)
        parser.add_argument("--gnomad-genomes-AF", type=float)

        parser.add_argument("--topmed-AC", type=int)
        parser.add_argument("--topmed-AF", type=float)

        # annotations
        parser.add_argument("--transcript-consequence", nargs="*", help="transcript consequence terms", choices=[])
        parser.add_argument("--clinvar", nargs="*", help="clinvar pathogenicity terms", choices=[d['slug'] for d in CLINVAR_ANNOTATION_DEFINITIONS])
        parser.add_argument("--clinvar-stars", help="minimum number of stars", choices=range(4))
        parser.add_argument("--hgmd", nargs="*", help="hgmd terms")

        # allele frequencies2
        parser.add_argument("--this-callset-AC2", type=int)
        parser.add_argument("--this-callset-AF2", type=float)

        parser.add_argument("--g1k-AC2", type=int)
        parser.add_argument("--g1k-AF2", type=float)

        parser.add_argument("--exac-AC2", type=int)
        parser.add_argument("--exac-HH2", type=int)
        parser.add_argument("--exac-AF2", type=float)

        parser.add_argument("--gnomad-exomes-AC2", type=int)
        parser.add_argument("--gnomad-exomes-HH2", type=int)
        parser.add_argument("--gnomad-exomes-AF2", type=float)

        parser.add_argument("--gnomad-genomes-AC2", type=int)
        parser.add_argument("--gnomad-genomes-HH2", type=int)
        parser.add_argument("--gnomad-genomes-AF2", type=float)

        parser.add_argument("--topmed-AC2", type=int)
        parser.add_argument("--topmed-AF2", type=float)

        # annotations2
        parser.add_argument("--transcript-consequence2", nargs="*", help="transcript consequence terms", choices=[])
        parser.add_argument("--clinvar2", nargs="*", help="clinvar pathogenicity terms", choices=[d['slug'] for d in CLINVAR_ANNOTATION_DEFINITIONS])
        parser.add_argument("--clinvar-stars2", help="minimum number of stars", choices=range(4))
        parser.add_argument("--hgmd2", nargs="*", help="hgmd terms")

    def handle(self, *args, **options):

        # parse args

        # select families
        project_category_guids = options['project-category']
        project_guids = options['projects']
        analysis_group_guids = options['analysis_group']
        family_guids = options['family']

        individual_guids = options['individual']  # limit search in this family to just these individuals

        # inheritance filters
        inheritance_mode = options["inheritance_mode"]
        allow_missing_genotypes = options["allow_missing_genotypes"]

        # QC filters
        qc_pass = options["pass"]
        ab_affected = options["AB_affected"]
        dp_affected = options["DP_affected"]
        gq_affected = options["GQ_affected"]

        # options["mq"]
        # options["qd"]

        # loci
        gene_lists = options["gene_list"]
        gene_ids = options["gene_id"]
        transcript_ids = options["transcript_id"]
        intervals = options["interval"]
        variants = options["variant"]

        #dataset type
        dataset_types = options["dataset_type"]
        sample_types = options["sample_type"]
        variant_types = options["variant_type"]

        # allele frequencies
        parser.add_argument("--this-callset-AC", type=int)
        parser.add_argument("--this-callset-AF", type=float)

        parser.add_argument("--g1k-AC", type=int)
        parser.add_argument("--g1k-AF", type=float)

        parser.add_argument("--exac-AC", type=int)
        parser.add_argument("--exac-HH", type=int)
        parser.add_argument("--exac-AF", type=float)

        parser.add_argument("--gnomad-exomes-AC", type=int)
        parser.add_argument("--gnomad-exomes-HH", type=int)
        parser.add_argument("--gnomad-exomes-AF", type=float)

        parser.add_argument("--gnomad-genomes-AC", type=int)
        parser.add_argument("--gnomad-genomes-HH", type=int)
        parser.add_argument("--gnomad-genomes-AF", type=float)

        parser.add_argument("--topmed-AC", type=int)
        parser.add_argument("--topmed-AF", type=float)


        # compute the set of samples to search
        query_filter = (
            Q(individual__family__project__projectcategory__guid__in=project_category_guids) |
            Q(individual__family__project__guid__in=project_guids) |
            Q(individual__family__guid__in=family_guids))

        query_filter &= Q(sample_status=Sample.SAMPLE_STATUS_LOADED)
        query_filter &= Q(elasticsearch_index__isnull=False)
        if sample_types:
            query_filter &= Q(sample_type__in=sample_types)

        samples = Sample.objects.filter(query_filter).select_related("project", 'individual', 'individual__set')

        if options['user']:
            # confirm the user has VIEW permissions for these families
            user = User.objects.get(user=options['user'])

            forbidden_families = families.filter(~Q(can_view_group__user=user))
            if forbidden_families:
                raise PermissionDenied("Access denied to {} families".format(len(forbidden_families)))

        # get all elasticsearch indices to search


        # use MultiSearch to create a separate query against each one



        # TODO add gene-level annotations: RVIS, NMD, constraint, pLI, mis_Z, omim


        #if projects_to_process:
        #    projects = Project.objects.filter(Q(name__in=projects_to_process) | Q(guid__in=projects_to_process))
        #    logging.info("Processing %s projects" % len(projects))
        #else:
        #    projects = Project.objects.filter(deprecated_project_id__isnull=False)
        #    logging.info("Processing all %s projects" % len(projects))


        #dataset_types2 = options["dataset_type2"]
        #sample_types2 = options["sample_type2"]
        #variant_types2 = options["variant_type2"]
