import logging
import os
import psutil
from deploy.utils.kubectl_utils import get_pod_name, run_in_pod, run
from deploy.utils.servctl_utils import check_kubernetes_context

logger = logging.getLogger(__name__)


def load_dataset(deployment_target, project_name, genome_version, sample_type, dataset_type, vcf, memory_to_use=None, cpu_limit=None, **kwargs):
    pod_name = get_pod_name('pipeline-runner', deployment_target=deployment_target)

    # run load command
    additional_load_command_args = "  ".join("--%s '%s'" % (key.lower().replace("_", "-"), value) for key, value in kwargs.items() if value is not None)

    run_locally = deployment_target == "minikube"
    if run_locally:
        vcf_name = os.path.basename(vcf)
        path_in_pod = "/data/{}".format(vcf_name)
        if os.path.isfile(vcf):
            run("kubectl cp '%(vcf)s' '%(pod_name)s:%(path_in_pod)s'" % locals()) # if local file path, copy file into pod
        elif vcf.startswith("http"):
            run_in_pod(pod_name, "wget -N %(vcf)s -O %(path_in_pod)s" % locals())
        vcf = path_in_pod

        total_memory = psutil.virtual_memory().total - 6*10**9  # leave 6Gb for other processes
        memory_to_use = "%sG" % (total_memory / 2 / 10**9) if memory_to_use is None else memory_to_use # divide available memory evenly between spark driver & executor
        cpu_limit = max(1, psutil.cpu_count() / 2) if cpu_limit is None else cpu_limit
        load_command = """/hail-elasticsearch-pipelines/run_hail_locally.sh \
            --driver-memory %(memory_to_use)s \
            --executor-memory %(memory_to_use)s \
            hail_scripts/v01/load_dataset_to_es.py \
                --cpu-limit %(cpu_limit)s \
                --genome-version %(genome_version)s \
                --project-guid %(project_name)s \
                --sample-type %(sample_type)s \
                --dataset-type %(dataset_type)s \
                --skip-validation \
                --exclude-hgmd \
                --vep-block-size 100 \
                --es-block-size 10 \
                --num-shards 1 \
                --max-samples-per-index 99 \
                %(additional_load_command_args)s \
                %(vcf)s
        """ % locals()

    else:
        load_command = """/hail-elasticsearch-pipelines/run_hail_on_dataproc.sh \
            hail_scripts/v01/load_dataset_to_es.py \
                --genome-version %(genome_version)s \
                --project-guid %(project_name)s \
                --sample-type %(sample_type)s \
                --dataset-type %(dataset_type)s \
                %(additional_load_command_args)s \
                %(vcf)s
        """ % locals()

    run_in_pod(pod_name, load_command, verbose=True)


def load_example_project(deployment_target, genome_version="37", cpu_limit=None, start_with_step=None):
    """Load example project

    Args:
        deployment_target (string):
        genome_version (string): reference genome version - either "37" or "38"
    """

    project_name = "1kg"

    check_kubernetes_context(deployment_target)

    pod_name = get_pod_name('seqr', deployment_target=deployment_target)
    if not pod_name:
        raise ValueError("No 'seqr' pod found. Is the kubectl environment configured in this terminal? and have either of these pods been deployed?" % locals())

    run_in_pod(pod_name, "wget -N https://storage.googleapis.com/seqr-reference-data/test-projects/1kg.ped" % locals())
    #run_in_pod(pod_name, "gsutil cp %(ped)s ." % locals())

    # TODO call APIs instead?
    run_in_pod(pod_name, "python2.7 -u -m manage create_project -p '1kg.ped' '%(project_name)s'" % locals(), verbose=True)

    if genome_version == "37":
        vcf_filename = "1kg.vep.vcf.gz"
    elif genome_version == "38":
        vcf_filename = "1kg.liftover.GRCh38.vep.vcf.gz"
    else:
        raise ValueError("Unexpected genome_version: %s" % (genome_version,))

    load_dataset(
        deployment_target,
        project_name=project_name,
        genome_version=genome_version,
        sample_type="WES",
        dataset_type="VARIANTS",
        cpu_limit=cpu_limit,
        start_with_step=start_with_step,
        vcf="https://storage.googleapis.com/seqr-reference-data/test-projects/%(vcf_filename)s" % locals())


def update_reference_data(deployment_target):
    """Load older reference data.

    Args:
        deployment_target (string):
    """

    check_kubernetes_context(deployment_target)

    pod_name = get_pod_name('seqr', deployment_target=deployment_target)
    if not pod_name:
        raise ValueError("No 'seqr' pods found. Is the kubectl environment configured in this terminal? and have either of these pods been deployed?" % locals())

    # commented out because this is not loaded from settings backup
    #run_in_pod(pod_name, "python2.7 -u manage.py update_all_reference_data --omim-key '$OMIM_KEY'" % locals(), verbose=True, print_command=True)

    run_in_pod(pod_name, "mkdir -p /seqr/data/reference_data")
    run_in_pod(pod_name, "wget -N https://storage.googleapis.com/seqr-reference-data/seqr-resource-bundle.tar.gz -O /seqr/data/reference_data/seqr-resource-bundle.tar.gz")
    run_in_pod(pod_name, "tar xzf /seqr/data/reference_data/seqr-resource-bundle.tar.gz -C /seqr/data/reference_data", verbose=True)
    run_in_pod(pod_name, "rm /seqr/data/reference_data/seqr-resource-bundle.tar.gz")

    # load legacy resources
    run_in_pod(pod_name, "python -u manage.py load_resources", verbose=True)
    run_in_pod(pod_name, "python -u manage.py load_omim", verbose=True)


def create_user(deployment_target, email=None, password=None):
    """Creates a seqr superuser

    Args:
        deployment_target (string):
        email (string): if provided, user will be created non-interactively
        password (string): if provided, user will be created non-interactively
    """
    check_kubernetes_context(deployment_target)

    if not email:
        run_in_pod("seqr", "python -u manage.py createsuperuser" % locals(), is_interactive=True)
    else:
        logger.info("Creating user %(email)s" % locals())
        run_in_pod("seqr",
            """echo "from django.contrib.auth.models import User; User.objects.create_superuser('%(email)s', '%(email)s', '%(password)s')" \| python manage.py shell""" % locals(),
            print_command=False,
            errors_to_ignore=["already exists"])
