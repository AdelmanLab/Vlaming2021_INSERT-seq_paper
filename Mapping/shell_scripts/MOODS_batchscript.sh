#!/bin/bash
#SBATCH -c 2                               # Request two cores
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-02:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=6000                         # Memory total in MB (for all cores)
#SBATCH -o hv34_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hv34_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=hanneke_vlaming@hms.harvard.edu   # Email to which notifications will be sent

module load python/2.7.12
cd ~/scratch3_hv34/MOODS/

source ~/adelman_folder/KW_working/MOODS/.envrcO2
python ~/adelman_folder/KW_working/MOODS/src/MOODS-python-1.9.3/scripts/moods_dna.py -m motifs/"$1".pfm -s ~/Library_design/Final_tables/sequence_inread_full.fa -p "$2" > output/MOODS_"$1"_output_p"$2"