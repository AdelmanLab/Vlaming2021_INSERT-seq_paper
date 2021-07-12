#!/bin/bash
#SBATCH -c 2                               # Request two cores
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-02:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4000                         # Memory total in MB (for all cores)
#SBATCH -o hv34_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hv34_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=hanneke_vlaming@hms.harvard.edu   # Email to which notifications will be sent

module load perl/5.30.0
module load homer/4.10.3

findMotifs.pl input/top_"$1"_topoverbottom_seqs.fa mouse output/"$1"_topoverbottom_RNA -fasta input/bottom_"$1"_topoverbottom_seqs.fa -norevopp -len 6,8,10 -mcheck input/ATtRACT.cb
