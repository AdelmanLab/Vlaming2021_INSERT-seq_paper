#!/bin/bash
#SBATCH -c 3                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-01:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=500                         # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=hanneke_vlaming@hms.harvard.edu   # Email to which notifications will be sent

# Assume I'm already inside /home/hv34/scratch3_hv34/mapping_sortseq/
# Assume I have 3 cores
module load bowtie2/2.3.4.3
module load samtools/1.9
	
# Assume I've already made the fake genome
#bowtie2-build --threads 3 ~/screens_analysis/Reference/STAR_October2020/sequences_unspliced_spliced_pooleddata.fa ~/screens_analysis/Reference/STAR_October2020/sequences_unspliced_spliced_pooleddata_fake_genome_bt2

# Assume I've already indexed the reference genome
#samtools faidx ~/screens_analysis/Reference/STAR_October2020/sequences_unspliced_spliced_pooleddata.fa

sbatch ~/scratch3_hv34/mapping_sortseq/Bowtie2_afterSTAR_onRNApersample.sh LIB042739_GEN00162228
sbatch ~/scratch3_hv34/mapping_sortseq/Bowtie2_afterSTAR_onRNApersample.sh LIB043895_GEN00166547
sbatch ~/scratch3_hv34/mapping_sortseq/Bowtie2_afterSTAR_onRNApersample.sh LIB043895_GEN00166548
sbatch ~/scratch3_hv34/mapping_sortseq/Bowtie2_afterSTAR_onRNApersample.sh LIB043895_GEN00166549
sbatch ~/scratch3_hv34/mapping_sortseq/Bowtie2_afterSTAR_onRNApersample.sh LIB043895_GEN00166550
sbatch ~/scratch3_hv34/mapping_sortseq/Bowtie2_afterSTAR_onRNApersample.sh t_10M

sbatch ~/scratch3_hv34/mapping_sortseq/Bowtie2_afterSTAR_onRNApersample.sh run_723
sbatch ~/scratch3_hv34/mapping_sortseq/Bowtie2_afterSTAR_onRNApersample.sh run_725