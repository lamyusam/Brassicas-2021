#!/bin/bash

# Request wallclock time (format hours:minutes:seconds).
#SBATCH --time=24:00:00

# Request RAM
# #SBATCH --mem=10G

# Request TMPDIR space.
# #SBATCH --tmp=10G

# Select the number of threads.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

# Set the name of the task.
#SBATCH --job-name=BLAST_raph_vs_arabid

# Set working directory
#SBATCH --chdir=/home/bat1m21/scratch/brassica-project/analysis/Raphanus/RBH_BLAST

#load BLAST
module load blast/2.11.0

# create BLAST databases
makeblastdb -in species_fastas/arabidopsis.faa -dbtype 'prot' -out databases/arabidopsis_db
makeblastdb -in species_fastas/raph.faa -dbtype 'prot' -out databases/raphanus_db

# run BLAST search
# -evalue: equivalent to significance value cutoff (lower values = more stringent matching)
# -outfmt: output format (6 = tabular)
blastp -num_threads 8 -evalue 1 -outfmt 6 -query species_fastas/arabidopsis.faa -db databases/raphanus_db > arabid_to_raph
blastp -num_threads 8 -evalue 1 -outfmt 6 -query species_fastas/raph.faa -db databases/arabidopsis_db > raph_to_arabid