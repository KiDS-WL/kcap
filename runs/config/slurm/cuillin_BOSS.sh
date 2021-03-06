#!/bin/bash
#SBATCH --ntasks 80
#SBATCH --mem-per-cpu 2000                         # In MB
#SBATCH --time 2-00:00                      # Time in days-hours:min
#SBATCH --job-name=kcap_BOSS               # this will be displayed if you write squeue in terminal and will be in the title of all emails slurm sends
#SBATCH --requeue                           # Allow requeing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ttr@roe.ac.uk
#SBATCH -o runs/output/BOSS/logs/log_fid.out
#SBATCH -e runs/output/BOSS/logs/log_fid.err
module load intel
module load openmpi
source ${HOME}/Codes/miniconda/bin/activate kcap_test

mpirun --host `python ${HOME}/Codes/SlurmEnvToHostfile/SlurmEnvToHostfile.py --no-file` --mca btl ^openib --oversubscribe cosmosis --mpi runs/config/BOSS_only.ini
