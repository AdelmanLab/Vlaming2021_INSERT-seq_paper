#!/bin/bash
#SBATCH -c 3                               # Request three cores
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-04:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5000                         # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=hanneke_vlaming@hms.harvard.edu   # Email to which notifications will be sent

module load star/2.7.0a

cd ./STAR_output_20210611_Oct2020/

STAR --runThreadN 3 --genomeDir ../Reference/STAR --readFilesIn ../cutadapt_output/cDNA_total2-4_R1_trim.cutadapt.fastq ../cutadapt_output/cDNA_total2-4_R2_trim.cutadapt.fastq --outSJfilterOverhangMin 10 4 4 4 --peOverlapNbasesMin 10 --peOverlapMMp 0.08 --alignEndsType EndToEnd --outFilterMismatchNmax 5 --scoreDelOpen -10 --scoreInsOpen -10
