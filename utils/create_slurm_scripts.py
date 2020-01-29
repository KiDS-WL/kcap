#!/usr/bin/env python
import os
import argparse
import sys
import pathlib

def create_sbatch_config(run_name, config_file, output_dir, log_dir, args, job_name=None, print_config=True):
    jobname = job_name or run_name
    template = f"""#!/bin/bash
#SBATCH --ntasks {args.n_task}
#SBATCH --mem-per-cpu 2000                         # In MB
#SBATCH --time {args.time}                      # Time in days-hours:min
#SBATCH --job-name={jobname}               # this will be displayed if you write squeue in terminal and will be in the title of all emails slurm sends
#SBATCH --requeue                           # Allow requeing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ttr@roe.ac.uk
#SBATCH -o {log_dir}/log_{run_name}.out
#SBATCH -e {log_dir}/log_{run_name}.err
#SBATCH --partition={args.partition}
module load intel
module load openmpi/4.0.0/intel
source ${{HOME}}/Codes/miniconda/bin/activate kcap_env

export OMP_NUM_THREADS={args.n_thread}

mpirun --host `python ${{HOME}}/Codes/SlurmEnvToHostfile/SlurmEnvToHostfile.py --no-file` --mca btl ^openib cosmosis --mpi {config_file}
    """
    if print_config: print(template)

    if not args.root_dir:
        with open(os.path.join(output_dir, "slurm_script_command.sh"), "w") as f:
            f.write(" ".join(sys.argv))
    with open(os.path.join(output_dir, "sbatch_command.sh"), "w") as f:
        f.write(template)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root-dir", required=True)

    parser.add_argument("--config_file")
    parser.add_argument("--run-name")
    parser.add_argument("--job-name")
    parser.add_argument("--output-dir")
    parser.add_argument("--log-dir", default="runs/logs")

    parser.add_argument("--n-task", default=32)
    parser.add_argument("--partition", default="all")
    parser.add_argument("--time", default="4-00:00")
    parser.add_argument("--n-thread", default="1")
    args = parser.parse_args()

    if args.root_dir:
        print(f"Generating sbatch files for all sub directories of root dir {args.root_dir}")
        path = pathlib.Path(args.root_dir)
        config_files = path.glob("*/config/pipeline.ini")
        if not config_files:
            print("No config files found.")

        for config_file in config_files:
            run_name = config_file.parts[-3]
            job_name = run_name
            output_dir = config_file.parent
            create_sbatch_config(run_name, config_file, output_dir, args.log_dir, args, job_name, print_config=False)
            print(f"Processed {run_name}.")

        sbatch_files = path.glob("*/config/sbatch_command.sh")
        with open(path.joinpath("sbatch_all.sh"), "w") as f:
            f.write("#!/bin/sh\n")
            f.writelines([f"sbatch {sbf}\n" for sbf in sbatch_files])
        
    else:
        create_sbatch_config(args.run_name, args.config_file, args.output_dir, args.log_dir, args, args.job_name)


