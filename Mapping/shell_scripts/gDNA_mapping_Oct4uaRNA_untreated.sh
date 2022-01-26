#!/bin/bash
#SBATCH -c 3                               # Request three cores
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-00:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=500                         # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=hanneke_vlaming@hms.harvard.edu   # Email to which notifications will be sent

# Assume I'm already inside /home/hv34/scratch3_hv34/mapping_sortseq/
# Assume I have 3 cores
module load bowtie2/2.3.4.3
module load cutadapt/1.14
module load samtools/1.9

# Assume I've already made and indexed the fake genome

sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB042739_GEN00162229
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB042739_GEN00162230
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB042739_GEN00162231
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB042739_GEN00162232
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB042739_GEN00162233
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB042739_GEN00162234
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB042739_GEN00162235
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB042739_GEN00162236
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB042739_GEN00162237

sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB043895_GEN00166541
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB043895_GEN00166542
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB043895_GEN00166543
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB043895_GEN00166544
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB043895_GEN00166545
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB043895_GEN00166546
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB043895_GEN00166551
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB043895_GEN00166552
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample.sh LIB043895_GEN00166553

sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample_novogene.sh gDNA_723
sbatch ~/scratch3_hv34/mapping_sortseq/gDNA_mapping_persample_novogene.sh gDNA_725