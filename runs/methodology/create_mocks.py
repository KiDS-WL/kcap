import sys
import subprocess

if __name__ == "__main__":

    if len(sys.argv) > 1:
        noise_begin = int(sys.argv[1])
        noise_end   = int(sys.argv[2])
    else:
        noise_begin = 0
        noise_end   = 1

    script = "utils/run_kcap.py"
    
    KiDS_twopoint_file =  "../K1000-data/Phase-1/twopoint/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad.fits"

    run_name_root = "base"
    # Create 3x2pt data vector
    run_type = "EE_nE_w"

    # Create noiseless mock data vector
    output_dir = "runs/methodology/data/noisefree_fiducial/"
    
    if noise_begin == 0:
        run_name = f"{run_name_root}_{run_type}"

        cmd = ["--create-mocks", "--noiseless-mocks",
                "--root-dir", output_dir,
                "--KiDS-data-file", KiDS_twopoint_file,
                "--run-name", run_name,
                "--run-type", run_type]
        subprocess.run(["python", script] + cmd, check=True)

    # Create noisy mocks
    output_dir = "runs/methodology/data/noisy_fiducial/"

    for i in range(noise_begin, noise_end):
        run_name = f"{run_name_root}_{i}_{run_type}"

        cmd = ["--create-mocks",
                "--root-dir", output_dir,
                "--KiDS-data-file", KiDS_twopoint_file,
                "--run-name", run_name,
                "--run-type", run_type]
        subprocess.run(["python", script] + cmd, check=True)


