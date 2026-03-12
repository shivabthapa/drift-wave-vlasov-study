#!/bin/bash
#SBATCH -J hw_turbulence
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH -o hw.out
#SBATCH -e hw.err

cd $SLURM_SUBMIT_DIR

/dartfs/rc/lab/E/EPaCO/libs/gkylsoft/gkeyll/bin/gkeyll drift_wave_instability.lua