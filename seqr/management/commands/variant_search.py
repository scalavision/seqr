import logging

import elasticsearch
import elasticsearch_dsl
import redis
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.core.management.base import BaseCommand
from django.db.models.query_utils import Q as Q_sql
from elasticsearch_dsl import Q as Q_es


from seqr.models import Project, Family, Sample, LocusListGene, LocusListInterval
from seqr.views.utils.permissions_utils import check_permissions
from xbrowse.core.constants import CLINVAR_ANNOTATION_DEFINITIONS

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = 'Variant search'

    def add_arguments(self, parser):

        parser.add_argument("--host", help="elasticsearch host", required=True)
        parser.add_argument("--port", help="elasticsearch port", default=9200)

        parser.add_argument("--genome-version", help="Search will be performed in this genome version", default=37)

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
        parser.add_argument("--dataset-type", nargs="*", choices=[Sample.DATASET_TYPE_VARIANT_CALLS, Sample.DATASET_TYPE_SV_CALLS])
        parser.add_argument("--sample-type", nargs="*", choices=[Sample.SAMPLE_TYPE_WES, Sample.SAMPLE_TYPE_WGS])
        parser.add_argument("--variant-type", nargs="*", choices=["SNP", "INS", "DEL"])

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
        options["this_callset_AC"]
        options["this_callset_AF"]

        options["g1k_AC"]
        options["g1k_AF"]

        options["exac_AC"]
        options["exac_HH"]
        options["exac_AF"]

        options["gnomad_exomes_AC"]
        options["gnomad_exomes_HH"]
        options["gnomad_exomes_AF"]

        options["gnomad_genomes_AC"]
        options["gnomad_genomes_HH"]
        options["gnomad_genomes_AF"]

        options["topmed_AC"]
        options["topmed_AF"]

        # compute the set of samples to search
        query_filter = (
            Q_sql(individual__family__project__projectcategory__guid__in=project_category_guids) |
            Q_sql(individual__family__project__guid__in=project_guids) |
            Q_sql(individual__family__guid__in=family_guids))

        query_filter &= Q_sql(sample_status=Sample.SAMPLE_STATUS_LOADED)
        query_filter &= Q_sql(elasticsearch_index__isnull=False)
        if sample_types:
            query_filter &= Q_sql(sample_type__in=sample_types)
        if dataset_types:
            query_filter &= Q_sql(dataset_type__in=dataset_types)

        samples = Sample.objects.filter(query_filter).select_related('individual', 'family', 'project')

        if options['user']:
            # confirm the user has VIEW permissions for these families
            user = User.objects.get(user=options['user'])

            forbidden_samples = samples.filter(~Q_sql(individual__family__project__can_view_group__user=user))
            if forbidden_samples:
                families = set(s.individual.family for s in forbidden_samples)
                projects = set(f.project.name for f in families)
                raise PermissionDenied("You don't have permissions to search the {} families in project(s): {}".format(
                    len(families), ", ".join(projects)))

        es_client = elasticsearch.Elasticsearch(host=options["host"])

        try:
            redis_client = redis.StrictRedis(host=options["redis-host"], socket_connect_timeout=3)
            redis_client.ping()
        except redis.exceptions.TimeoutError as e:
            logger.warn("Unable to connect to redis host: {}".format(options["redis-host"]) + str(e))
            redis_client = None

        # build site-level ES filters



        #~Q('term', **{encoded_sample_id+"_num_alt": 1}) |
        #Q('range', **{encoded_sample_id+"_ab": {'gte': min_ab}}))
        #s = s.filter(~Q('exists', field='filters'))
        #s = s.filter('range', **{encoded_sample_id+"_gq": {'gte': min_gq}})
        # build inheritance query for each family

        # get all elasticsearch indices




        # use MultiSearch to create a separate query against each one
        multi_search = elasticsearch_dsl.MultiSearch(using=es_client)


        # TODO add gene-level annotations: RVIS, NMD, constraint, pLI, mis_Z, omim


def _build_elasticsearch_site_filters(gene_list_guids=(), pass_only=False):
    site_filter = Q_es('term')

    # QC
    if pass_only:
        site_filter &= ~Q_es('exists', field='filters')

    if gene_list_guids:
        genes = LocusListGene.objects.filter(locus_list__guid__in=gene_list_guids).only('gene_id')
        site_filter &= Q_es('terms', [g.gene_id for g in genes])

        # intervals = options["interval"]
        intervals = LocusListInterval.objects.filter(locus_list__guid__in=gene_list_guids).only('genome_version', 'chrom', 'start', 'end')

        # run liftover and cache results
        

    # loci
    gene_lists = options["gene_list"]
    gene_ids = options["gene_id"]
    transcript_ids = options["transcript_id"]


    variants = options["variant"]


"""
elasticsearch index fields:

###  chrom-pos-ref-alt
contig  	string
pos  	number					
rsid  	string					
xpos  	number					
xstart  	number					
xstop  	number
start  	number					
end  	number					
ref  	string					
alt  	string
originalAltAlleles  	string					
wasSplit  	boolean										
aIndex  	number					
docId  	string					
variantId  	string		

####  site-level QC
filters  	string					
InbreedingCoeff  	number					
MQ  	number					
QD  	number				

###  misc 				
_id  	string
_index  	string					
_score  	number					
_source  	_source					
_type  	string					


###  VEP
geneIds  	string					
codingGeneIds  	string
transcriptConsequenceTerms  	string					
transcriptIds  	string
domains  	string		
		
sortedTranscriptConsequences  	string					

mainTranscript_amino_acids  	string					
mainTranscript_biotype  	string					
mainTranscript_canonical  	number					
mainTranscript_category  	string					
mainTranscript_cdna_end  	number					
mainTranscript_cdna_start  	number					
mainTranscript_codons  	string					
mainTranscript_domains  	string					
mainTranscript_gene_id  	string					
mainTranscript_gene_symbol  	string					
mainTranscript_hgvs  	string					
mainTranscript_hgvsc  	string					
mainTranscript_hgvsp  	string					
mainTranscript_lof  	string					
mainTranscript_lof_filter  	string					
mainTranscript_lof_flags  	string					
mainTranscript_major_consequence  	string					
mainTranscript_major_consequence_rank  	number					
mainTranscript_protein_id  	string					
mainTranscript_transcript_id  	string


###  clinvar / HGMD
clinvar_allele_id  	number					
clinvar_clinical_significance  	string					
clinvar_gold_stars  	number					
clinvar_review_status  	string
hgmd_accession  	string					
hgmd_class  	string					


###  deleteriousness predictors		
cadd_PHRED  	number						
dbnsfp_DANN_score  	number					
dbnsfp_FATHMM_pred  	string					
dbnsfp_GERP_RS  	string					
dbnsfp_MetaSVM_pred  	string					
dbnsfp_MutationTaster_pred  	string					
dbnsfp_Polyphen2_HVAR_pred  	string					
dbnsfp_REVEL_score  	string					
dbnsfp_SIFT_pred  	string					
dbnsfp_phastCons100way_vertebrate  	string								
eigen_Eigen_phred  	number					
primate_ai_score  	number					
mpc_MPC  	string	



###  allele frequencies
AC  	number					
AF  	number					
AN  	number

exac_AC_Adj  	number					
exac_AC_Hemi  	number					
exac_AC_Het  	number					
exac_AC_Hom  	number					
exac_AF  	number					
exac_AF_POPMAX  	number					
exac_AN_Adj  	number					

g1k_AC  	number					
g1k_AF  	number					
g1k_AN  	number					
g1k_POPMAX_AF  	number					

gnomad_exome_coverage  	number					
gnomad_exomes_AC  	number					
gnomad_exomes_AF  	number					
gnomad_exomes_AF_POPMAX_OR_GLOBAL  	number					
gnomad_exomes_AN  	number					
gnomad_exomes_Hemi  	number					
gnomad_exomes_Hom  	number				
	
gnomad_genome_coverage  	number					
gnomad_genomes_AC  	number					
gnomad_genomes_AF  	number					
gnomad_genomes_AF_POPMAX_OR_GLOBAL  	number					
gnomad_genomes_AN  	number					
gnomad_genomes_Hemi  	number					
gnomad_genomes_Hom  	number					

topmed_AC  	number					
topmed_AF  	number					
topmed_AN  	number					
topmed_Het  	number					
topmed_Hom  	number	


###  genotypes
sample0_ab  	number					
sample0_dp  	number					
sample0_gq  	number					
sample0_num_alt  	number					

sample1_ab  	number					
sample1_dp  	number					
sample1_gq  	number					
sample1_num_alt  	number					

"""