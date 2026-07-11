#!/bin/bash

#SBATCH --time 3-00:00

#SBATCH --output="logs/snakemake/%A-%x-std-output.out"
#SBATCH --error="logs/snakemake/%A-%x-err-output.out"

#SBATCH --mail-user=daiki.tagami@hertford.ox.ac.uk
#SBATCH --mail-type=ALL

#SBATCH --mem=10G
#SBATCH --cpus-per-task=1

#SBATCH --job-name="snakemake-slurm"

#SBATCH --cluster=swan
#SBATCH --partition=standard-statgen-cpu
#SBATCH --nodelist=smew01.cpu.stats.ox.ac.uk

snakemake --workflow-profile profiles --use-conda
