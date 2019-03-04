#!/bin/bash
#SBATCH -J all_anvio_proka
#SBATCH --mem=234G                   # Memory requested in megabytes. If omitted, the default is 1024 MB.
#SBATCH -t 7-00:00:00      # How long will your job run for? If omitted, the default is 3 hours. days-hours:minutes:seconds
#SBATCH -c 24

#trim adapters, and quality control with Trimmomatic
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_16.GR_DNA_1_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_16.GR_DNA_1_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_16.GR_DNA_1.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_16.GR_DNA_1_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_16.GR_DNA_1_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_16.GR_DNA_1.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_16.GR_DNA_1_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_16.GR_DNA_1_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_16.GR_DNA_1.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_28.GR_DNA_3_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_28.GR_DNA_3_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_28.GR_DNA_3.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_28.GR_DNA_3_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_28.GR_DNA_3_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_28.GR_DNA_3.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_28.GR_DNA_3_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_28.GR_DNA_3_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_28.GR_DNA_3.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_40.GR_DNA_5_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_40.GR_DNA_5_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_40.GR_DNA_5.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_40.GR_DNA_5_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_40.GR_DNA_5_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_40.GR_DNA_5.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_40.GR_DNA_5_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_40.GR_DNA_5_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_40.GR_DNA_5.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_52.GR_DNA_11_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_52.GR_DNA_11_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_52.GR_DNA_11.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_52.GR_DNA_11_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_52.GR_DNA_11_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_52.GR_DNA_11.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_52.GR_DNA_11_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_52.GR_DNA_11_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_52.GR_DNA_11.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_64.GR_DNA_14_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_64.GR_DNA_14_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_64.GR_DNA_14.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_64.GR_DNA_14_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_64.GR_DNA_14_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_64.GR_DNA_14.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_64.GR_DNA_14_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_64.GR_DNA_14_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_64.GR_DNA_14.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_76.GR_DNA_17_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_76.GR_DNA_17_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.001.BioOHT_76.GR_DNA_17.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_76.GR_DNA_17_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_76.GR_DNA_17_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.002.BioOHT_76.GR_DNA_17.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_76.GR_DNA_17_R1.fastq.gz ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_76.GR_DNA_17_R2.fastq.gz -baseout ~/matti/Hazen_metagenome/HI.4261.003.BioOHT_76.GR_DNA_17.fastq.gz ILLUMINACLIP:/global/home/hpc3229/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:3:30:8:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:30 AVGQUAL:20

#concatenate reads for the sample
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_16.GR_DNA_1_1P.fastq.gz > ~/matti/Hazen_metagenome/DNA_1_1P.fastq.gz
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_16.GR_DNA_1_2P.fastq.gz > ~/matti/Hazen_metagenome/DNA_1_2P.fastq.gz
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_16.GR_DNA_3_1P.fastq.gz > ~/matti/Hazen_metagenome/DNA_3_1P.fastq.gz
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_16.GR_DNA_3_2P.fastq.gz > ~/matti/Hazen_metagenome/DNA_3_2P.fastq.gz
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_16.GR_DNA_5_1P.fastq.gz > ~/matti/Hazen_metagenome/DNA_5_1P.fastq.gz
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_16.GR_DNA_5_2P.fastq.gz > ~/matti/Hazen_metagenome/DNA_5_2P.fastq.gz
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_16.GR_DNA_11_1P.fastq.gz > ~/matti/Hazen_metagenome/DNA_11_1P.fastq.gz
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_16.GR_DNA_11_2P.fastq.gz > ~/matti/Hazen_metagenome/DNA_11_2P.fastq.gz
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_16.GR_DNA_14_1P.fastq.gz > ~/matti/Hazen_metagenome/DNA_14_1P.fastq.gz
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_16.GR_DNA_14_2P.fastq.gz > ~/matti/Hazen_metagenome/DNA_14_2P.fastq.gz
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_16.GR_DNA_17_1P.fastq.gz > ~/matti/Hazen_metagenome/DNA_17_1P.fastq.gz
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_16.GR_DNA_17_2P.fastq.gz > ~/matti/Hazen_metagenome/DNA_17_2P.fastq.gz

#concatenate forward and reverse reads for matam
cat ~/matti/Hazen_metagenome/DNA_1_*P.fastq.gz > ~/matti/Hazen_metagenome/DNA_1.fastq.gz
cat ~/matti/Hazen_metagenome/DNA_3_*P.fastq.gz > ~/matti/Hazen_metagenome/DNA_3.fastq.gz
cat ~/matti/Hazen_metagenome/DNA_5_*P.fastq.gz > ~/matti/Hazen_metagenome/DNA_5.fastq.gz
cat ~/matti/Hazen_metagenome/DNA_11_*P.fastq.gz > ~/matti/Hazen_metagenome/DNA_11.fastq.gz
cat ~/matti/Hazen_metagenome/DNA_14_*P.fastq.gz > ~/matti/Hazen_metagenome/DNA_14.fastq.gz
cat ~/matti/Hazen_metagenome/DNA_17_*P.fastq.gz > ~/matti/Hazen_metagenome/DNA_17.fastq.gz
gunzip ~/matti/Hazen_metagenome/DNA_1.fastq.gz
gunzip ~/matti/Hazen_metagenome/DNA_3.fastq.gz
gunzip ~/matti/Hazen_metagenome/DNA_5.fastq.gz
gunzip ~/matti/Hazen_metagenome/DNA_11.fastq.gz
gunzip ~/matti/Hazen_metagenome/DNA_14.fastq.gz
gunzip ~/matti/Hazen_metagenome/DNA_17.fastq.gz

#run matam to find and classify full-length 16S sequences
source activate py36
matam_assembly.py -d ~/matti/SILVA/SILVA_128_SSURef_NR95 -i ~/matti/Hazen_metagenome/DNA_1.fastq -o ~/matti/Hazen_metagenome/matam/1/ --cpu 24 --max_memory 233000 -v --perform_taxonomic_assignment
matam_assembly.py -d ~/matti/SILVA/SILVA_128_SSURef_NR95 -i ~/matti/Hazen_metagenome/DNA_3.fastq -o ~/matti/Hazen_metagenome/matam/3/ --cpu 24 --max_memory 233000 -v --perform_taxonomic_assignment
matam_assembly.py -d ~/matti/SILVA/SILVA_128_SSURef_NR95 -i ~/matti/Hazen_metagenome/DNA_5.fastq -o ~/matti/Hazen_metagenome/matam/5/ --cpu 24 --max_memory 233000 -v --perform_taxonomic_assignment
matam_assembly.py -d ~/matti/SILVA/SILVA_128_SSURef_NR95 -i ~/matti/Hazen_metagenome/DNA_11.fastq -o ~/matti/Hazen_metagenome/matam/11/ --cpu 24 --max_memory 233000 -v --perform_taxonomic_assignment
matam_assembly.py -d ~/matti/SILVA/SILVA_128_SSURef_NR95 -i ~/matti/Hazen_metagenome/DNA_14.fastq -o ~/matti/Hazen_metagenome/matam/14/ --cpu 24 --max_memory 233000 -v --perform_taxonomic_assignment
matam_assembly.py -d ~/matti/SILVA/SILVA_128_SSURef_NR95 -i ~/matti/Hazen_metagenome/DNA_17.fastq -o ~/matti/Hazen_metagenome/matam/17/ --cpu 24 --max_memory 233000 -v --perform_taxonomic_assignment

#concatenate reads for all the samples
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_*.GR_DNA_*_1P.fastq.gz > ~/matti/Hazen_metagenome/DNA_all_1P.fastq.gz
cat ~/matti/Hazen_metagenome/HI.4261.00*.BioOHT_*.GR_DNA_*_2P.fastq.gz > ~/matti/Hazen_metagenome/DNA_all_2P.fastq.gz

#quality check with seqtk and fastqc
seqtk fqchk ~/matti/Hazen_metagenome/DNA_all_1P.fastq.gz > ~/matti/Hazen_metagenome/fastqc/DNA_all_1P_seqtk.txt
seqtk fqchk ~/matti/Hazen_metagenome/DNA_all_2P.fastq.gz > ~/matti/Hazen_metagenome/fastqc/DNA_all_2P_seqtk.txt
fastqc ~/matti/Hazen_metagenome/DNA_all_1P.fastq.gz -o ~/matti/Hazen_metagenome/fastqc/
fastqc ~/matti/Hazen_metagenome/DNA_all_2P.fastq.gz -o ~/matti/Hazen_metagenome/fastqc/

#run Megahit
source activate py36
megahit -1 ~/matti/Hazen_metagenome/DNA_all_1P.fastq.gz -2 ~/matti/Hazen_metagenome/DNA_all_2P.fastq.gz -o ~/matti/Hazen_metagenome/megahit/all -t 24

#prepare the prokarya db
source activate anvio4
anvi-script-reformat-fasta ~/matti/Hazen_metagenome/megahit/all/final.contigs.fa -o ~/matti/Hazen_metagenome/megahit/all/final_contigs_anvio.fa -l 1000 --simplify-names
bowtie2-build ~/matti/Hazen_metagenome/megahit/all/final_contigs_anvio.fa ~/matti/Hazen_metagenome/mapping/all_proka/contigs --threads 24

#read coverage
bowtie2 --threads 24 -x ~/matti/Hazen_metagenome/mapping/all_proka/contigs -1 ~/matti/Hazen_metagenome/DNA_1_1P.fastq.gz -2 ~/matti/Hazen_metagenome/DNA_1_2P.fastq.gz --no-unal -S ~/matti/Hazen_metagenome/mapping/DNA_1_proka.sam
bowtie2 --threads 24 -x ~/matti/Hazen_metagenome/mapping/all_proka/contigs -1 ~/matti/Hazen_metagenome/DNA_3_1P.fastq.gz -2 ~/matti/Hazen_metagenome/DNA_3_2P.fastq.gz --no-unal -S ~/matti/Hazen_metagenome/mapping/DNA_3_proka.sam
bowtie2 --threads 24 -x ~/matti/Hazen_metagenome/mapping/all_proka/contigs -1 ~/matti/Hazen_metagenome/DNA_5_1P.fastq.gz -2 ~/matti/Hazen_metagenome/DNA_5_2P.fastq.gz --no-unal -S ~/matti/Hazen_metagenome/mapping/DNA_5_proka.sam
bowtie2 --threads 24 -x ~/matti/Hazen_metagenome/mapping/all_proka/contigs -1 ~/matti/Hazen_metagenome/DNA_11_1P.fastq.gz -2 ~/matti/Hazen_metagenome/DNA_11_2P.fastq.gz --no-unal -S ~/matti/Hazen_metagenome/mapping/DNA_11_proka.sam
owtie2 --threads 24 -x ~/matti/Hazen_metagenome/mapping/all_proka/contigs -1 ~/matti/Hazen_metagenome/DNA_14_1P.fastq.gz -2 ~/matti/Hazen_metagenome/DNA_14_2P.fastq.gz --no-unal -S ~/matti/Hazen_metagenome/mapping/DNA_14_proka.sam
bowtie2 --threads 24 -x ~/matti/Hazen_metagenome/mapping/all_proka/contigs -1 ~/matti/Hazen_metagenome/DNA_17_1P.fastq.gz -2 ~/matti/Hazen_metagenome/DNA_17_2P.fastq.gz --no-unal -S ~/matti/Hazen_metagenome/mapping/DNA_17_proka.sam
samtools view -F 4 -bS ~/matti/Hazen_metagenome/mapping/DNA_1_proka.sam > ~/matti/Hazen_metagenome/mapping/DNA_1_proka-RAW.bam
samtools view -F 4 -bS ~/matti/Hazen_metagenome/mapping/DNA_3_proka.sam > ~/matti/Hazen_metagenome/mapping/DNA_3_proka-RAW.bam
samtools view -F 4 -bS ~/matti/Hazen_metagenome/mapping/DNA_5_proka.sam > ~/matti/Hazen_metagenome/mapping/DNA_5_proka-RAW.bam
samtools view -F 4 -bS ~/matti/Hazen_metagenome/mapping/DNA_11_proka.sam > ~/matti/Hazen_metagenome/mapping/DNA_11_proka-RAW.bam
samtools view -F 4 -bS ~/matti/Hazen_metagenome/mapping/DNA_14_proka.sam > ~/matti/Hazen_metagenome/mapping/DNA_14_proka-RAW.bam
samtools view -F 4 -bS ~/matti/Hazen_metagenome/mapping/DNA_17_proka.sam > ~/matti/Hazen_metagenome/mapping/DNA_17_proka-RAW.bam
anvi-init-bam ~/matti/Hazen_metagenome/mapping/DNA_1_proka-RAW.bam -o ~/matti/Hazen_metagenome/mapping/DNA_1_proka.bam
anvi-init-bam ~/matti/Hazen_metagenome/mapping/DNA_3_proka-RAW.bam -o ~/matti/Hazen_metagenome/mapping/DNA_3_proka.bam
anvi-init-bam ~/matti/Hazen_metagenome/mapping/DNA_5_proka-RAW.bam -o ~/matti/Hazen_metagenome/mapping/DNA_5_proka.bam
anvi-init-bam ~/matti/Hazen_metagenome/mapping/DNA_11_proka-RAW.bam -o ~/matti/Hazen_metagenome/mapping/DNA_11_proka.bam
anvi-init-bam ~/matti/Hazen_metagenome/mapping/DNA_14_proka-RAW.bam -o ~/matti/Hazen_metagenome/mapping/DNA_14_proka.bam
anvi-init-bam ~/matti/Hazen_metagenome/mapping/DNA_17_proka-RAW.bam -o ~/matti/Hazen_metagenome/mapping/DNA_17_proka.bam
rm ~/matti/Hazen_metagenome/mapping/DNA_*_proka.sam ~/matti/Hazen_metagenome/mapping/DNA*_proka-RAW.bam ~/matti/Hazen_metagenome/mapping/*_proka.bt2

#calculate average read lengths for each sample for normalization
cat ~/matti/Hazen_metagenome/DNA_1_1P.fastq.gz ~/matti/Hazen_metagenome/DNA_1_2P.fastq.gz > ~/matti/Hazen_metagenome/fastqc/DNA_1_seqtk.fastq.gz
cat ~/matti/Hazen_metagenome/DNA_3_1P.fastq.gz ~/matti/Hazen_metagenome/DNA_3_2P.fastq.gz > ~/matti/Hazen_metagenome/fastqc/DNA_3_seqtk.fastq.gz
cat ~/matti/Hazen_metagenome/DNA_5_1P.fastq.gz ~/matti/Hazen_metagenome/DNA_5_2P.fastq.gz > ~/matti/Hazen_metagenome/fastqc/DNA_5_seqtk.fastq.gz
cat ~/matti/Hazen_metagenome/DNA_11_1P.fastq.gz ~/matti/Hazen_metagenome/DNA_11_2P.fastq.gz > ~/matti/Hazen_metagenome/fastqc/DNA_11_seqtk.fastq.gz
cat ~/matti/Hazen_metagenome/DNA_14_1P.fastq.gz ~/matti/Hazen_metagenome/DNA_14_2P.fastq.gz > ~/matti/Hazen_metagenome/fastqc/DNA_14_seqtk.fastq.gz
cat ~/matti/Hazen_metagenome/DNA_17_1P.fastq.gz ~/matti/Hazen_metagenome/DNA_17_2P.fastq.gz > ~/matti/Hazen_metagenome/fastqc/DNA_17_seqtk.fastq.gz
seqtk fqchk ~/matti/Hazen_metagenome/fastqc/DNA_1_seqtk.fastq.gz | sed -n "1p" | awk -F ";" -v OFS='\t' '{print 1, $3}' > ~/matti/Hazen_metagenome/fastqc/DNA_1_seqtk.tmp
seqtk fqchk ~/matti/Hazen_metagenome/fastqc/DNA_3_seqtk.fastq.gz | sed -n "1p" | awk -F ";" -v OFS='\t' '{print 3, $3}' > ~/matti/Hazen_metagenome/fastqc/DNA_3_seqtk.tmp
seqtk fqchk ~/matti/Hazen_metagenome/fastqc/DNA_5_seqtk.fastq.gz | sed -n "1p" | awk -F ";" -v OFS='\t' '{print 5, $3}' > ~/matti/Hazen_metagenome/fastqc/DNA_5_seqtk.tmp
seqtk fqchk ~/matti/Hazen_metagenome/fastqc/DNA_11_seqtk.fastq.gz | sed -n "1p" | awk -F ";" -v OFS='\t' '{print 11, $3}' > ~/matti/Hazen_metagenome/fastqc/DNA_11_seqtk.tmp
seqtk fqchk ~/matti/Hazen_metagenome/fastqc/DNA_14_seqtk.fastq.gz | sed -n "1p" | awk -F ";" -v OFS='\t' '{print 14, $3}' > ~/matti/Hazen_metagenome/fastqc/DNA_14_seqtk.tmp
seqtk fqchk ~/matti/Hazen_metagenome/fastqc/DNA_17_seqtk.fastq.gz | sed -n "1p" | awk -F ";" -v OFS='\t' '{print 17, $3}' > ~/matti/Hazen_metagenome/fastqc/DNA_17_seqtk.tmp
cat ~/matti/Hazen_metagenome/fastqc/DNA_*_seqtk.tmp > ~/matti/Hazen_metagenome/fastqc/sample_average_read_lengths.txt
rm ~/matti/Hazen_metagenome/fastqc/*seqtk.*

#formulate the database
anvi-gen-contigs-database -f ~/matti/Hazen_metagenome/megahit/all/final_contigs_anvio.fa -o ~/matti/Hazen_metagenome/anvio/all_proka/contigs.db -n 'LakeHazen prokaryotes'

#continue the regular anvio pipeline
hmm profiles
anvi-run-hmms -c ~/matti/Hazen_metagenome/anvio/all_proka/contigs.db --num-threads 24
#NCBI cogs
anvi-run-ncbi-cogs -c ~/matti/Hazen_metagenome/anvio/all_proka/contigs.db --num-threads 24
#create profiles for each sample
for sample in `awk 'NR > 1 {print $1}' ~/matti/Hazen_metagenome/samples_proka.txt`; do anvi-profile -i ~/matti/Hazen_metagenome/mapping/$sample.bam -c ~/matti/Hazen_metagenome/anvio/all_proka/contigs.db -o ~/matti/Hazen_metagenome/anvio/$sample --num-threads 24 --min-contig-length 1000; done
#merge everything and bin genomes
anvi-merge ~/matti/Hazen_metagenome/anvio/DNA*/PROFILE.db -o ~/matti/Hazen_metagenome/anvio/SAMPLES-MERGED -c ~/matti/Hazen_metagenome/anvio/all_proka/contigs.db --sample-name 'LakeHazen'

#manual refinement of bins with anvi-refine

#annotate the anvio database CDS (after splitting with http://iubio.bio.indiana.edu/gmod/genogrid/scripts/split_multifasta.pl) for function with InterProScan
anvi-get-aa-sequences-for-gene-calls -c ~/matti/Hazen_metagenome/anvio/refined/all_proka/contigs.db -o ~/matti/Hazen_metagenome/anvio/refined/all_proka/protein-sequences.fa
mkdir -p ~/matti/Hazen_metagenome/anvio/refined/all_proka/split_sequences/
split_multifasta.pl -i ~/matti/Hazen_metagenome/anvio/refined/all_proka/protein-sequences.fa -o ~/matti/Hazen_metagenome/anvio/refined/all_proka/split_sequences/ --seqs_per_file=5000

#run an array of 372 protein multifastas as individual profiles in IPRS
sbatch ~/matti/Hazen_metagenome/IPRS_array_proka.sh

#format all the outputs for anvi-import-functions so that it works (export COGs, and re-import with the others). Tool used to parse the GO annotations: https://github.com/xvazquezc/stuff/raw/master/iprsalt2anvio.sh
cat ~/matti/Hazen_metagenome/anvio/refined/all_proka/split_sequences/interpro-output_*.tsv > ~/matti/Hazen_metagenome/anvio/refined/all_proka/interpro-output.tsv
cut -f1,4,5,6,9 ~/matti/Hazen_metagenome/anvio/refined/all_proka/interpro-output.tsv > ~/matti/Hazen_metagenome/anvio/refined/all_proka/interpro-output_anvio.tsv
iprsalt2anvio.sh ~/matti/Hazen_metagenome/anvio/refined/all_proka/interpro-output.tsv ~/matti/Hazen_metagenome/anvio/refined/all_proka/interpro-output_anvio_GO.tsv
head -n -1 ~/matti/Hazen_metagenome/anvio/refined/all_proka/interpro-output_anvio_GO.tsv | tail -n +2 > temp.tsv ; mv temp.tsv ~/matti/Hazen_metagenome/anvio/refined/all_proka/interpro-output_anvio_GO.tsv
cat ~/matti/Hazen_metagenome/anvio/refined/all_proka/split_sequences/*_hits_table.txt | grep -v "#" | awk -v OFS="\t" '{print $3, "Resfams", $2, $1, $5}' > ~/matti/Hazen_metagenome/anvio/refined/all_proka/resfams_annotations.tsv
anvi-export-functions -c ~/matti/Hazen_metagenome/anvio/refined/all_proka/contigs.db -o ~/matti/Hazen_metagenome/anvio/refined/all_proka/COG_annotations.tsv --annotation-sources COG_FUNCTION,COG_CATEGORY
cat ~/matti/Hazen_metagenome/anvio/refined/all_proka/COG_annotations.tsv ~/matti/Hazen_metagenome/anvio/refined/all_proka/interpro-output_anvio_GO.tsv ~/matti/Hazen_metagenome/anvio/refined/all_proka/interpro-output_anvio.tsv ~/matti/Hazen_metagenome/anvio/refined/all_proka/resfams_annotations.tsv > temp.tsv
mv temp.tsv ~/matti/Hazen_metagenome/anvio/refined/all_proka/interpro-output_anvio.tsv
anvi-import-functions -c ~/matti/Hazen_metagenome/anvio/refined/all_proka/contigs.db -i ~/matti/Hazen_metagenome/anvio/refined/all_proka/interpro-output_anvio.tsv --drop-previous-annotations

#summarize the collection
anvi-summarize -p ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-MERGED/PROFILE.db -c ~/matti/Hazen_metagenome/anvio/refined/all_proka/contigs.db -C CONCOCT -o ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ --init-gene-coverages

#check and summarize the quality and completion of the bins with CheckM
mkdir -p ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/all_bins/
cp -s ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/*/*-contigs.fa ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/all_bins/
source activate qiime1
mkdir -p ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/checkm/
checkm lineage_wf -x fa -t 24 ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/all_bins/ ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/checkm/
checkm qa -t 24 ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/checkm/lineage.ms ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/checkm/ -f ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/checkm/qa.out
cat ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/checkm/qa.out | sed 's/^[ \t]*//;s/[ \t]*$//' > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/checkm/qa_trimmed.out

#gather lineage info and make a genome tree of the bins
checkm tree_qa ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/checkm/ -o 2 -f ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/checkm/tree_qa_lineage.out --tab_table
checkm tree_qa ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/checkm/ -o 3 -f ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/checkm/tree_qa.tre

#calculate MetaCyc and KEGG pathways with MinPath and antismash in an array
sbatch ~/matti/Hazen_metagenome/Metagenome_array_proka.sh

#combine files for analysis
cat ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT/bin_by_bin/*/*kegg_pathways.txt > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT/kegg_pathways.txt
cat ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT/bin_by_bin/*/*metacyc_pathways.txt > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT/metacyc_pathways.txt
cat ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT/bin_by_bin/*/*antismash_clusters.txt > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT/antismash_clusters.txt
tail -q -n +2 ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT/bin_by_bin/*/*gene_coverages.txt > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT/gene_coverages.txt

#gather a table of individual marker genes
for annotation in `awk '{print $1}' ~/matti/Hazen_metagenome/individual_gene_annotations.txt`; do grep $annotation ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT/bin_by_bin/Bin_*/*-gene_calls.txt | cut -f1,2,3,4,5 | awk -v x=${annotation} -F "/|\t|:" '{print x, $11, $13, $14, $15, $16, $17}' > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT/individual_marker_genes/${annotation}_gene_calls.tmp & done; wait
cat ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT/individual_marker_genes/*_gene_calls.tmp > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT/individual_marker_genes/all_gene_calls.txt

#fetch the ribosomal proteins for the bins, align them with translatorX and concatenate alignments
mkdir -p ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/
for bin in $(awk '{print $0}' ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/quality_bins.txt)
do for annotation in `awk '{print $1}' ~/matti/Hazen_metagenome/RP_annotations.txt`
do grep $annotation ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/${bin}/${bin}-gene_calls.txt | awk -v OFS="" -v x=${annotation} -v y=${bin} '{print ">", x, "_", y, "\n", $NF}' > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${bin}_${annotation}.fa &
done
wait
done

#align the RPs with translatorx
for annotation in $(awk '{print $1}' ~/matti/Hazen_metagenome/RP_annotations.txt)
do
(cat ~/ncbi/GbBac/ribosomal_proteins/${annotation}.fa ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/*_${annotation}.fa > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}.fa
perl ~/.local/bin/translatorx_vLocal_custom.pl -i ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}.fa -o  ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_mafft_translatorx -c 11 -p F) &
done
wait

#remove the outliers, realign and trim with trimal
for annotation in $(awk '{print $1}' ~/matti/Hazen_metagenome/RP_annotations.txt)
do
(faSomeRecords ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_mafft_translatorx.nt_ali.fasta ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_bridging_sequences.list ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_bridging_sequences.fasta
faSomeRecords -exclude ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_mafft_translatorx.nt_ali.fasta ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_bridging_sequences.list ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_refined_mafft_translatorx.nt_ali.fasta
sed -i 's/-//g' ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_refined_mafft_translatorx.nt_ali.fasta
perl ~/.local/bin/translatorx_vLocal_custom.pl -i ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_refined_mafft_translatorx.nt_ali.fasta -o ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_refined_v2_mafft_translatorx -c 11 -p F
~/.local/bin/trimal -in ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_refined_v2_mafft_translatorx.nt_ali.fasta -out ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal.fasta -automated1

#remove duplicates (more than one copy of RP in a genome), correct headers and construct approximate ML trees with FastTree
for annotation in $(awk '{print $1}' ~/matti/Hazen_metagenome/RP_annotations.txt)
do
(fasuniq -i ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal.fasta | awk -F " " '/^>/{$0=$1}1' > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal_unique.fasta
~/.local/bin/FastTreeMP -gtr -gamma -nt ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal_unique.fasta  > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal.tre
treeshrink.py -q 0.01 -i ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal.tre
tr -s "\t" "\n" < ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal_treeshrink/${annotation}_trimal_shrunk_RS_0.01.txt > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal_treeshrink/${annotation}_outliers.list
faSomeRecords -exclude ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal_unique.fasta ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal_treeshrink/${annotation}_outliers.list ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal_final.fasta
~/.local/bin/FastTreeMP -gtr -gamma -nt ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal_final.fasta > ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/${annotation}_trimal_final.tre) &
done
wait

#extract 16S sequences for each bin
for bin in $(tail -n +2 ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bins_summary.txt | awk '{print $1}')
do
(ssu-align -f --no-align ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/${bin}/${bin}-contigs.fa ~/matti/Hazen_metagenome/anvio/refined/SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/${bin}/ssu-align/ ) &
done
wait