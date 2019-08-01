#!/usr/bin/env python

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("--job-name", default="kcap")
    parser.add_argument("--n-task", default=32)
    parser.add_argument("--partition", default="all")

    args = parser.parse_args()

    if "KV450" in args.config and "BOSS" in args.config:
        outputdir = "KV450_BOSS"
    elif "KV450" in args.config:
        outputdir = "KV450"
    elif "BOSS" in args.config:
        outputdir = "BOSS"
    else:
        raise ValueError("This doesn't seem to be either KV450 nor BOSS. Aborting.")

    template = f"""#!/bin/bash
#SBATCH --ntasks {args.n_task}
#SBATCH --mem-per-cpu 2000                         # In MB
#SBATCH --time 2-00:00                      # Time in days-hours:min
#SBATCH --job-name={args.job_name}               # this will be displayed if you write squeue in terminal and will be in the title of all emails slurm sends
#SBATCH --requeue                           # Allow requeing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ttr@roe.ac.uk
#SBATCH -o runs/output/{outputdir}/logs/log_{args.job_name}.out
#SBATCH -e runs/output/{outputdir}/logs/log_{args.job_name}.err
#SBATCH --partition={args.partition}
module load intel
module load openmpi
source ${{HOME}}/Codes/miniconda/bin/activate kcap_test

mpirun --host `python ${{HOME}}/Codes/SlurmEnvToHostfile/SlurmEnvToHostfile.py --no-file` --mca btl ^openib --oversubscribe cosmosis --mpi {args.config}
    """
    print(template)
