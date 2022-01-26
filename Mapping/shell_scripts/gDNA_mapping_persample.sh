#!/bin/bash
#SBATCH -c 3                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-04:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4000                         # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=hanneke_vlaming@hms.harvard.edu   # Email to which notifications will be sent

# Assume I've already made the fake genome
# Assume I've already indexed the reference genome
# Assume I've already loaded modules and have defined variable 1

# Trim adapters from read, and have quality cutoff
cutadapt -f fastq --match-read-wildcards -m 20 -q 20,20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ~/scratch3_hv34/mapping_sortseq/cutadapt_output/"$1"_R1_trim.cutadapt.fastq -p ~/scratch3_hv34/mapping_sortseq/cutadapt_output/"$1"_R2_trim.cutadapt.fastq ~/scratch3_hv34/mapping_sortseq/fastq_files/"$1"_R1.fastq ~/scratch3_hv34/mapping_sortseq/fastq_files/"$1"_R2.fastq > ~/scratch3_hv34/mapping_sortseq/cutadapt_output/"$1"_cutadapt_log.txt

# Run allowing ~3.2 mismatches (so: max 3 at perfect bases, maybe 4 of imperfect bases), disallow gaps completely
bowtie2 --end-to-end -p 3 --no-discordant --no-mixed --rdg 20,5 --rfg 20,5 --score-min L,0,-0.11 --un ~/scratch3_hv34/mapping_sortseq/mapped_count_files/"$1"_bt2_PE_unhit_3wogap -x ~/scratch3_hv34/mapping_sortseq/Reference/inserts_fake_genome_bt2 -1 ~/scratch3_hv34/mapping_sortseq/cutadapt_output/"$1"_R1_trim.cutadapt.fastq -2 ~/scratch3_hv34/mapping_sortseq/cutadapt_output/"$1"_R2_trim.cutadapt.fastq -S ~/scratch3_hv34/mapping_sortseq/mapped_count_files/"$1"_bt2_PE_nodisc_3mm_wogaps

# Make counts table from the bowtie2 result
samtools import ~/scratch3_hv34/mapping_sortseq/Reference/sequence_inread_full.fa.fai ~/scratch3_hv34/mapping_sortseq/mapped_count_files/"$1"_bt2_PE_nodisc_3mm_wogaps ~/scratch3_hv34/mapping_sortseq/mapped_count_files/"$1"_bt2_PE_nodisc_3mm_wogaps.bam
samtools sort ~/scratch3_hv34/mapping_sortseq/mapped_count_files/"$1"_bt2_PE_nodisc_3mm_wogaps.bam -o ~/scratch3_hv34/mapping_sortseq/mapped_count_files/"$1"_bt2_PE_nodisc_3mm_wogaps_sorted.bam
samtools index ~/scratch3_hv34/mapping_sortseq/mapped_count_files/"$1"_bt2_PE_nodisc_3mm_wogaps_sorted.bam
samtools idxstats ~/scratch3_hv34/mapping_sortseq/mapped_count_files/"$1"_bt2_PE_nodisc_3mm_wogaps_sorted.bam > ~/scratch3_hv34/mapping_sortseq/mapped_count_files/"$1"_bt2_PE_nodisc_3mm_wogaps_idxstats.txt