#!/usr/bin/env python
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--job-name")
    parser.add_argument("--n-task", default=32)
    parser.add_argument("--partition", default="all")
    parser.add_argument("--time", default="2-00:00")
    parser.add_argument("--n-thread", default="1")
    args = parser.parse_args()

    outputdir = os.path.join(args.run_dir, args.run_name)
    jobname = args.job_name or args.run_name
    template = f"""#!/bin/bash
#SBATCH --ntasks {args.n_task}
#SBATCH --mem-per-cpu 2000                         # In MB
#SBATCH --time {args.time}                      # Time in days-hours:min
#SBATCH --job-name={jobname}               # this will be displayed if you write squeue in terminal and will be in the title of all emails slurm sends
#SBATCH --requeue                           # Allow requeing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ttr@roe.ac.uk
#SBATCH -o {outputdir}/logs/log_{args.run_name}.out
#SBATCH -e {outputdir}/logs/log_{args.run_name}.err
#SBATCH --partition={args.partition}
module load intel
module load openmpi/4.0.0/intel
source ${{HOME}}/Codes/miniconda/bin/activate kcap_env

export OMP_NUM_THREADS={args.n_thread}

mpirun --host `python ${{HOME}}/Codes/SlurmEnvToHostfile/SlurmEnvToHostfile.py --no-file` --mca btl ^openib cosmosis --mpi {args.config}
    """
    print(template)
