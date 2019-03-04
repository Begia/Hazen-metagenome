#!/bin/bash
#SBATCH -J IPRS_array
#SBATCH --mem=16G       # Memory requested in megabytes. If omitted, the default is 1024 MB.
#SBATCH -t 1-00:00:00      # How long will your job run for? If omitted, the default is 3 hours. days-hours:minutes:seconds
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=1-372
#SBATCH --output=IPRSarrayJob_%A_%a.out

source activate anvio4
interproscan.sh -i ~/matti/Hazen_metagenome/anvio/refined/all_proka/split_sequences/${SLURM_ARRAY_TASK_ID}.fsa -f tsv -o ~/matti/Hazen_metagenome/anvio/refined/all_proka/split_sequences/interpro-output_${SLURM_ARRAY_TASK_ID}.tsv -goterms -pa -appl Pfam,TIGRFAM,PIRSF