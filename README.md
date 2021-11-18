# Hazen metagenome
Code used in the analysis of data in the article, "Microbial genomes retrieved from High Arctic lake sediments encode for adaptation to cold and oligotrophic environments" by Matti O. Ruuskanen, Graham Colby, Kyra A. St.Pierre, Vincent L. St.Louis,Stéphane Aris-Brosou & Alexandre J. Poulain (https://doi.org/10.1002/lno.11334).

## Included files

1. Main sequence handling shell script
	- Scripts that were run to produce the primary data from Illumina read files
	- One file script: Primary_shell_scripts.sh

2. Accessory shell scripts (used by the sequence handling script)
	- IPRS_array_proka.sh
		* Runs an array of 372 protein multifastas as individual profiles in interproscan in parallel
	- Metagenome_array_proka.sh
		* Calculates MetaCyc and KEGG pathways with MinPath and antismash in an array (all MAGs in parallel)
		
3. Data analysis scripts
	- Data_analysis.R
	- Scripts that use the primary data to produce all the figures, tables, and other results in the manuscript
