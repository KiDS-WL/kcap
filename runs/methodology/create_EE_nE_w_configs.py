import subprocess
import os

if __name__ == "__main__":

    script = "utils/run_kcap.py"
        
    # Need to change?
    dz_cov_file = "data/KV450/nofz/DIR_cov.asc"

    # Case for 3x2pt. For the full setup we'd want to loop over the run_type
    run_type  = "EE_nE_w"

    # MAP runs
    output_root_dir = "runs/methodology/chains/MAP"
    run_name_root = "MAP"

    n_noise_mocks = 1

    for i in range(n_noise_mocks):
        root_data_dir = f"runs/methodology/data/noisy_fiducial/base_{i}_EE_nE_w/data/"
        twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgretta_nOfZ_bucerosBroad_mock_noisy.fits")
        boss_data_files = [os.path.join(root_data_dir, "BOSS/BOSS_mock_noisy_bin_1.txt"),
                           os.path.join(root_data_dir, "BOSS/BOSS_mock_noisy_bin_2.txt")]
        boss_cov_files =  [os.path.join(root_data_dir, "BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt"),
                           os.path.join(root_data_dir, "BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt")]

        run_name = f"{run_name_root}_{i}_{run_type}" 
        cmd = ["--root-dir", output_root_dir,
               "--run-name", run_name,
               "--run-type", run_type,
               "--KiDS-data-file", twopoint_file,
               "--dz-covariance-file", dz_cov_file,
               "--BOSS-data-files", *boss_data_files,
               "--BOSS-covariance-files", *boss_cov_files,
               "--sampler", "maxlike"]
        subprocess.run(["python", script] + cmd, check=True)

    # Multinest and test sampler runs
    root_data_dir = "runs/methodology/data/noisefree_fiducial/base_EE_nE_w/data/"
    twopoint_file = os.path.join(root_data_dir, "KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgretta_nOfZ_bucerosBroad_mock_noiseless.fits")
    boss_data_files = [os.path.join(root_data_dir, "BOSS/BOSS_mock_noiseless_bin_1.txt"),
                       os.path.join(root_data_dir, "BOSS/BOSS_mock_noiseless_bin_2.txt")]
    boss_cov_files =  [os.path.join(root_data_dir, "BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt"),
                       os.path.join(root_data_dir, "BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt")]

    # Multinest
    output_root_dir = "runs/methodology/chains/multinest"
    multinest_settings = ["--sampler-config", "multinest_efficiency", "0.3",
                          "--sampler-config", "nested_sampling_tolerance", "1.0e-2"]

    run_name_root = "mulitnest"
    run_name = f"{run_name_root}_{run_type}" 
    cmd = ["--root-dir", output_root_dir,
            "--run-name", run_name,
            "--run-type", run_type,
            "--KiDS-data-file", twopoint_file,
            "--dz-covariance-file", dz_cov_file,
            "--BOSS-data-files", *boss_data_files,
            "--BOSS-covariance-files", *boss_cov_files,
            "--sampler", "multinest",
            *multinest_settings]
    subprocess.run(["python", script] + cmd, check=True)

    # Test sampler
    output_root_dir = "runs/methodology/chains/test"
    run_name_root = "test_sampler"
    run_name = f"{run_name_root}_{run_type}" 
    cmd = ["--root-dir", output_root_dir,
            "--run-name", run_name,
            "--run-type", run_type,
            "--KiDS-data-file", twopoint_file,
            "--dz-covariance-file", dz_cov_file,
            "--BOSS-data-files", *boss_data_files,
            "--BOSS-covariance-files", *boss_cov_files,
            "--sampler", "test"]
    subprocess.run(["python", script] + cmd, check=True)
