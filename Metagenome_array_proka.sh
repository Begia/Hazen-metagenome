#!/bin/bash
#SBATCH -J MetagenomeArray
#SBATCH --mem=8G       # Memory requested in megabytes. If omitted, the default is 1024 MB.
#SBATCH -t 00-08:00:00      # How long will your job run for? If omitted, the default is 3 hours. days-hours:minutes:seconds
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=2-1489
#SBATCH --output=MetagenomeArrayJob_%A_%a.out

BIN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bins_summary.txt | cut -f1)

source activate qiime1
grep -o 'GO:[0-9]\{7\}' ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-gene_calls.txt | sort | uniq  > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-terms.go
awk 'FNR==NR{go2ec[$NF]=$1; next} { if ($0 in go2ec) { ec=go2ec[$0]; n=split(ec, t, "."); ecshow=substr(ec, 4); for(i=n+1;i<=4;i++) ecshow=ecshow".-"; print $1 " " ecshow; }}' ~/bin/MinPath/data/ec2go ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-terms.go > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-terms.ec
MinPath1.4.py -any ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-terms.ec -map ~/bin/MinPath/data/ec2path -report ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-metacyc.minpath -details ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-metacyc.details
awk -v x=$BIN -v OFS='\t' 'FNR==NR{if ($8 == 1) pwy[$14]=$14; next} ($1 in pwy) { print x, $0}' ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-metacyc.minpath ~/bin/MinPath/data/pathway.info > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-metacyc_pathways.txt

MinPath1.4.py -any ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-terms.ec -map ~/bin/MinPath/data/ec2kegg_final -report ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-kegg.minpath -details ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-kegg.details
awk -v x=$BIN -v OFS='\t' 'FNR==NR{if ($8 == 1) pwy[$14]=$14; next} ($1 in pwy) {print x, $0}' ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-kegg.minpath ~/bin/MinPath/data/kegg_pathway.info > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/$BIN/$BIN-kegg_pathways.txt