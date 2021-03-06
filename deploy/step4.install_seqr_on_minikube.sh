#!/usr/bin/env bash

echo ==== Wait for minikube to start =====
set -x

for i in {1..150}; do    # timeout for 5 minutes
   kubectl get po &> /dev/null
   if [ $? -ne 1 ]; then
      break
  fi
  echo 'Waiting for minikube to start...'
  sleep 2
done

sudo minikube addons disable dashboard  # disable kuberentes dashboard to conserve resources

minikube status


set +x
echo ==== deploy all seqr components =====

source venv/bin/activate

wget -N https://storage.googleapis.com/seqr-reference-data/gene_reference_data_backup.gz
./servctl deploy-all --restore-seqr-db-from-backup gene_reference_data_backup.gz --docker-image-tag release_20180926 minikube
