
seqr
====
[![Build Status](https://travis-ci.org/macarthur-lab/seqr.svg?branch=master)](https://travis-ci.org/macarthur-lab/seqr)

seqr is a web-based analysis tool for rare disease genomics.

This repository contains code that underlies the [Broad seqr instance](http://seqr.broadinstitute.org) and other seqr deployments.

## Technical Overview

seqr consists of the following components:
- seqr - the main client-server application. It consists of javascript + react.js on the client-side, python + django on the server-side.
- postgres - SQL database used by seqr and phenotips to store project metadata and user-generated content such as variant notes, etc.
- phenotips - 3rd-party web-based form for entering structured phenotype data.
- matchbox - a service that encapsulates communication with the Match Maker Exchange.
- nginx - http server used as the main gateway between seqr and the internet.
- pipeline-runner - container for running hail pipelines to annotate and load new datasets.
- redis - in-memory cache used to speed up request handling.
- elasticsearch - NoSQL database used to store variant callsets.
- kibana - dashboard and visual interface for elasticsearch.
- mongo - legacy NoSQL database originally used for variant callsets and still used now to store some reference data and logs.

These components can be deployed on local or cloud-based hardware.
A laptop with 4 CPUs and 16G RAM may be sufficient for small datasets.
The seqr [production instance](http://seqr.broadinstitute.org) currently uses two n1-highmem-4 (4 vCPUs, 26 GB memory) servers on google cloud + separate servers for the elasticsearch database.
The pipeline for loading new datasets uses Spark to parallelelize VEP and other annotation steps. This pipeline can run on the same machine that's hosting seqr components, but a separate Spark cluster will make the loading process faster.


## Install

Whether installing seqr on a laptop, on-prem, or on cloud VMs, we use Docker images and Kubernetes to automate the deployment steps and isolate them from the host systems.
Users that are very familiar with components that make up seqr may want to install them directly on host systems without using Kubernetes. This provides maximum control, but requires more work for installation and maintenance. 
If you do decide to go this route, the [Dockerfiles](https://github.com/macarthur-lab/seqr/tree/master/deploy/docker) can be useful as a list of steps for installing each component.  
  In most cases, we recommend the docker/kubernetes approach because:  1) it automates deployment as much as possible 2) as seqr evolves, it makes it easier to roll out new components or move around existing components if they outgrow existing hardware.   

The instructions below cover local deployments using Minikube, but are also directly applicable to cloud-based deployments which replace Minikube with a managed Kubernetes cluster like Google Container Engine.
The Kubernetes project maintains a [complete list of deployment options](https://kubernetes.io/docs/setup/pick-right-solution/).

#### Prerequisites
 - *Hardware:*  At least **16 Gb RAM**, **2 CPUs**, **50 Gb disk space**  

 - *Software:*  
     - python2.7    
     - on MacOS only: [homebrew](http://brew.sh/) package manager  
     - on Linux only: root access with sudo (because we run minikube without a virtualization layer, using --vm-driver=none)
    
#### Step 1: Install Kubernetes

Local and on-prem installations can use [MiniKube](https://kubernetes.io/docs/setup/minikube/) to create a self-contained kubernetes cluster on a single machine. 

Run the following command to install `gcc`, `java1.8`, `minikube`, `kubectl` and their dependencies on the host machine (using `brew`, `yum` or `apt-get`):

###### MacOS
```
SCRIPT=step1.macos.install_dependencies.sh && curl -L http://raw.githubusercontent.com/macarthur-lab/seqr/master/deploy/$SCRIPT -o $SCRIPT && chmod 777 $SCRIPT && source $SCRIPT
```

###### CentOS7 / RedHat

```
SCRIPT=step1.linux-centos7.install_dependencies.sh && curl -L http://raw.githubusercontent.com/macarthur-lab/seqr/master/deploy/$SCRIPT -o $SCRIPT && chmod 777 $SCRIPT && source $SCRIPT
```

###### Ubuntu 

```
SCRIPT=step1.linux-ubuntu18.install_dependencies.sh && curl -L http://raw.githubusercontent.com/macarthur-lab/seqr/master/deploy/$SCRIPT -o $SCRIPT && chmod 777 $SCRIPT && source $SCRIPT
```

#### Step 2: Install and start elasticsearch

Run this command to start an elasticsearch instance in the current directory: 
```
SCRIPT=step2.install_and_start_elasticsearch.sh && curl -L http://raw.githubusercontent.com/macarthur-lab/seqr/master/deploy/$SCRIPT -o $SCRIPT && chmod 777 $SCRIPT && source $SCRIPT
```

#### Step 3: Download seqr deployment scripts

In a new terminal, run this command to download seqr deployment scripts:
```
SCRIPT=step3.download_deployment_scripts.sh && curl -L http://raw.githubusercontent.com/macarthur-lab/seqr/master/deploy/$SCRIPT -o $SCRIPT && chmod 777 $SCRIPT && source $SCRIPT
```
  
Optionally edit deployment settings before proceeding to step 4:

* *deploy/kubernetes/minikube-settings.yaml* - contains settings like $MINIKUBE_DISK_SIZE.
* *deploy/secrets/shared/** - directories that contain keys, passwords and other sensitive info that shouldn't be shared publicly. You may at some point want to edit:
  * *gcloud/service-account-key.json* - allows gcloud and kubectl to access google cloud resources in your project from within pods. We provide a placeholder key which can access public resources.
  * *nginx/tls.cert* and *nginx/tls.key* - ssh keys that allow https access and avoid web browser "insecure website" warnings. https connections are critical for encrypting seqr logins, so you will want to order your own keys before making your seqr instance visible over the internet.
  * *seqr/postmark_server_token* - seqr uses this to send outgoing emails via postmark.com mail service
  * *seqr/omim_key* - api key for downloading the latest omim files. We provide a placeholder key, but you'll want to use your own.
  * *matchbox/nodes.json* - contains the list of all nodes that matchbox can connect to on the MME network, along with the authentication token for each node.


#### Step 4: Install seqr on minikube

Run this command to deploy all seqr components to minikube and load reference data:

```
SCRIPT=step4.install_seqr_on_minikube.sh && curl -L http://raw.githubusercontent.com/macarthur-lab/seqr/master/deploy/$SCRIPT -o $SCRIPT && chmod 777 $SCRIPT && source $SCRIPT
```

This step may run for 1 to 1.5 hours depending on internet speed. 


#### Create superuser 

After this finishes, you can create a super-user account by running these commands in the `seqr` subdirectory:

```
source ./activate_virtualenv.sh
  
./servctl create-user minikube 
```

#### Open seqr in browser

Then, to open seqr in your browser, you can do: 

```
open http://$(minikube ip):30003   
```

assuming you are on MacOS and minikube is running on the same machine as the web browser.

If minikube is running on a remote machine and you want to access seqr over the internet, you will need to first make sure 
traffic is being forwarded from the remote machine's http port (port 80) to seqr's Node port in minikube (30003). One way to do this is to start 
an ssh tunnel on the machine that's running minikube:

```
sudo ssh -v -i ~/.ssh/id_rsa -N -L 0.0.0.0:80:localhost:30003 ${USER}@$(hostname)
```   

NOTE: If you encounter `Permission denied (publickey)` errors, you may need to [generate ssh keys](https://help.github.com/enterprise/2.14/user/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/#generating-a-new-ssh-key) and do `cat ~/.ssh/id_rsa.pub >> ~/.ssh/authorized_keys` before starting the tunnel.   


## Load example project

To create an example seqr project and load an exome dataset with 16 individuals from the 1000 genomes project, run: 

```
source ./activate_virtualenv.sh

./servctl load-example-project --cpu-limit 1 minikube
```

This can take 16 to 24 hours on a machine with 16 Gb RAM and using 1 CPU for annotation. 
If more memory is available, raising the `--cpu-limit` will allow multiple instances of VEP to run in parallel (with each VEP instance needing ~4 Gb additional RAM), 
and parallelize all other steps in the pipeline, with a proportional reduction in runtime.     


## Creating projects and loading datasets

A seqr project groups together users that are collaborating on the analysis of one or more datasets. It encapsulates the variant data, pedigree information, and any tags or notes that users create during analysis.
  
To create a new project:  
1. login to seqr and click on `+ Create Project` in the bottom right.  
2. click on the new project. Also, note the project GUID in the url `/project/${project_guid}/project_page` - for example `R0003_demo_project1`
3. click `Edit Families and Individuals` and use the Bulk Upload form to upload a pedigree file in [.fam format](https://www.cog-genomics.org/plink2/formats#fam)  
   
To annotate and load a new dataset, run the `servctl load-dataset` command. For example: 
```
source ./activate_virtualenv.sh

./servctl load-dataset minikube --genome-version 37 --project-guid "${project_guid}" --sample-type WES --dataset-type VARIANTS --cpu-limit 1 --input-vcf ${vcf_path} 
``` 

where `${vcf_path}` is replaced with `/local/path/to/your_data.vcf.gz` or `gs://some-google-storage-bucket/path/to/your_data.vcf.gz`, and `${project_guid}` comes from the project page url.

Once the dataset finishes loading, you can add it to a seqr project using the `Edit Datasets` form on the Project page. 


## Updating / Migrating an older xBrowse Instance

[Update/Migration Instructions](https://github.com/macarthur-lab/seqr/blob/master/deploy/MIGRATE.md)


## Deploying and managing seqr components


The `./servctl` wrapper script provides sub-commands for deploying and interacting with seqr components running on kubernetes.

Run `./servctl -h` to see all available subcommands. The most commonly used ones are:

```
  ./servctl

      deploy-all  {deployment-target}                       # end-to-end deployment - deploys all seqr components and loads reference data
      deploy {component-name} {deployment-target}           # deploy some specific component(s)
          --restore-seqr-db-from-backup  <seqr_db_backup_file.sql.gz>   # deploying seqr with this option will also load the given seqrdb backup into postgres
          --restore-phenotips-db-from-backup  <xwiki_db_backup_file.sql.gz>   # deploying seqr with this option will also load the given phenotips xwikidb backup into postgres
      status {deployment-target}                            # print status of all kubernetes and docker subsystems
      set-env {deployment-target}                           # deploy one or more components
      dashboard {deployment-target}                         # open the kubernetes dasbhoard in a browser

      shell {component-name} {deployment-target}            # open a bash shell inside one of the component containers
      logs {component-name} {deployment-target}             # show logs for one or more components
      troubleshoot {component-name} {deployment-target}     # print more detailed info that may be useful for discovering why a component is failing during pod initialization
      connect-to {component-name} {deployment-target}       # shows logs, and also sets up a proxy so that the server running inside this component can be accessed from http://localhost:<port>

      copy-to {component-name} {deployment-target} {local-path}           # copy a local file to one of the pods
      copy-from {component-name} {deployment-target} {path} {local-path}  # copy a file from one of the pods to a local directory

      delete {component-name} {deployment-target}           # undeploys the component

      create-user {deployment-target}                       # create seqr superuser
      update-reference-data {deployment-target}             # update gene-level reference data in seqr
      load-example-project {deployment-target}              # run commands in seqr and pipeline-runner to load an example project


  *** {deployment-target} is one of:  minikube, gcloud-dev, or gcloud-prod
  *** {component-name} is one of: init-cluster, settings, secrets, nginx, postgres, phenotips, seqr,
                                  redis, pipeline-runner, matchbox, cockpit, kibana,
                                  external-mongo-connector, external-elasticsearch-connector,
                                  mongo, elasticsearch,
                                  es-client, es-master, es-data, es-kibana
```


## Data loading pipelines

seqr uses [hail](http://hail.is)-based pipelines to run VEP and add in other reference data before loading them into elasticsearch.
These pipelines can be run locally on a single machine or on-prem spark cluster, or on a cloud-based spark cluster like Google Dataproc.
We are working on integrating these pipelines so that they are launched and managed by seqr.
For now, they must be run manually, as shown in the examples below. The code for these pipelines is in [Data annotation and loading pipelines](https://github.com/macarthur-lab/hail-elasticsearch-pipelines)
and is automatically installed in the `pipeline-runner` component which is deployed as part of standard seqr deployment.

Example with seqr deployed to google cloud GKE, and using Google Dataproc to run the pipeline:
```
# these commands should be run locally on your laptop
git clone git@github.com:macarthur-lab/hail-elasticsearch-pipelines.git

cd hail-elasticsearch-pipelines
export IP=56.4.0.3   # IP address of elasticsearch instance running on google cloud on the project-default network so that it's visible to dataproc nodes, but not to the public internet.
export SEQR_PROJECT_GUID=R003_seqr_project3  # guid of existing seqr project
export GS_VCF_PATH=gs://seqr-datasets/GRCh38/my-new-dataset.vcf.gz   # VCF path on cloud storage

# this will create a new dataproc cluster and submit the pipeline to it
python gcloud_dataproc/load_GRCh38_dataset.py --host $IP --project-guid SEQR_PROJECT_GUID --sample-type WGS --dataset-type VARIANTS $GS_VCF_PATH --es-block-size 50

# after the pipeline completes successfully, you can link the new elasticsearch index to the seqr project by using the 'Edit Datasets' dialog on the project page.
```


## Kubernetes Resources


- Official Kuberentes User Guide:  https://kubernetes.io/docs/user-guide/
- 15 Kubernetes Features in 15 Minutes: https://www.youtube.com/watch?v=o85VR90RGNQ
- Kubernetes: Up and Running: https://www.safaribooksonline.com/library/view/kubernetes-up-and/9781491935668/
- The Children's Illustrated Guide to Kubernetes: https://deis.com/blog/2016/kubernetes-illustrated-guide/
