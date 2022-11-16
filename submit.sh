#!/bin/bash

#SBATCH -o MVM_24Threads_YES10%ReadErrorRate_10power2VirusHostRatio_Simulation1.log-%j
#SBATCH -n 75

source /etc/profile # Initialize the module command first

module load anaconda/2020a # Load Anaconda and MPI module


# MVM 3.08 x 10^2 virus to host ratio
python3 -B ont_virus_detection_simulation.py MVM 'Minute virus of mice.fasta' 50 'GCF_000223135.1_CriGri_1.0_genomic.fna' 1 '/home/raeuf/raeuf_notebook/bioinformatics_tools/apps/blast/db/RefSeq' 10power2 50 500 yes 10 24 3



