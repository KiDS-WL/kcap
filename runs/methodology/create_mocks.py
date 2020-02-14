import subprocess

if __name__ == "__main__":

    script = "utils/run_kcap.py"
    
    KiDS_twopoint_file =  "../K1000-data/Phase-1/twopoint/twoPoint_PneE+PeeE_mean_None_cov_theoryEgretta_nOfZ_bucerosBroad.fits"

    run_name_root = "base"
    # Create 3x2pt data vector
    run_type = "EE_nE_w"

    # Create noiseless mock data vector
    output_dir = "runs/methodology/data/noisefree_fiducial/"
    run_name = f"{run_name_root}_{run_type}"

    cmd = ["--create-mocks", "--noiseless-mocks",
            "--root-dir", output_dir,
            "--KiDS-data-file", KiDS_twopoint_file,
            "--run-name", run_name,
            "--run-type", run_type]
    subprocess.run(["python", script] + cmd, check=True)

    # Create noisy mocks
    output_dir = "runs/methodology/data/noisy_fiducial/"

    n_noise_mocks = 1

    for i in range(n_noise_mocks):
        run_name = f"{run_name_root}_{i}_{run_type}"

        cmd = ["--create-mocks",
                "--root-dir", output_dir,
                "--KiDS-data-file", KiDS_twopoint_file,
                "--run-name", run_name,
                "--run-type", run_type]
        subprocess.run(["python", script] + cmd, check=True)


