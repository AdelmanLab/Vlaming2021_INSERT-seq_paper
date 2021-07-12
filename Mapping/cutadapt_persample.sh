#!/bin/bash
#SBATCH -c 2                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-00:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=3000                         # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=hanneke_vlaming@hms.harvard.edu   # Email to which notifications will be sent


# Trim adapters from read, and have quality cutoff
module load cutadapt/1.14
cutadapt -f fastq --match-read-wildcards -m 20 -q 20,20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o /home/hv34/scratch3_hv34/mapping_sortseq/cutadapt_output/"$1"_R1_trim.cutadapt.fastq -p /home/hv34/scratch3_hv34/mapping_sortseq/cutadapt_output/"$1"_R2_trim.cutadapt.fastq /home/hv34/scratch3_hv34/mapping_sortseq/fastq_files/"$1"_R1.fastq /home/hv34/scratch3_hv34/mapping_sortseq/fastq_files/"$1"_R2.fastq > /home/hv34/scratch3_hv34/mapping_sortseq/cutadapt_output/"$1"_cutadapt_log.txt
