#!/bin/bash
#SBATCH --job-name=drift64
#SBATCH --partition=preemptable
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=drift64.out
#SBATCH --error=drift64.err

echo "Starting job at $(date)"

# Path to Gkeyll
GKYL=/dartfs/rc/lab/E/EPaCO/libs/gkylsoft/gkeyll/bin/gkeyll

srun $GKYL drift_wave_instability.lua

echo "Finished at $(date)"