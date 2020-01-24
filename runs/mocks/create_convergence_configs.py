import subprocess

if __name__ == "__main__":

    # utils/run_kcap.py --create-mocks --noiseless-mocks --root-dir data/KiDS1000/mocks/noisefree/ --KiDS-data-file data/KiDS1000/mocks/twoPoint_PneE+PeeE_mean_theoryNoiseFree_cov_theoryBuceros_nOfZ_bucerosBroad.fits --run-type EE_nE_w --run-name EE_nE_w_base --overwrite

    script = "utils/run_kcap.py"
    root_dir = "runs/mocks/convergence_tests"
    
    base_twopoint_file = "data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/KiDS/twoPoint_PneE+PeeE_mean_theoryNoiseFree_cov_theoryBuceros_nOfZ_bucerosBroad_mock_noiseless.fits"
    
    # Need to change
    base_dz_cov_file = "data/KV450/nofz/DIR_cov.asc"

    base_boss_data_files = ["data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS_mock_noiseless_bin_1.txt",
                            "data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS_mock_noiseless_bin_2.txt"]
    
    # base_boss_cov_files = ["data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt",
    #                        "data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt"]

    # Main chains

    twopoint_file = base_twopoint_file
    dz_cov_file = base_dz_cov_file
    boss_data_files = base_boss_data_files
    run_type = "EE_nE_w"

    sampler = "emcee"
    run_name_root = "emcee"
    run_name = f"{run_name_root}_{run_type}"
    cmd = ["--root-dir", root_dir,
           "--run-name", run_name,
           "--run-type", run_type,
           "--KiDS-data-file", twopoint_file,
           "--dz-covariance-file", dz_cov_file,
           "--BOSS-data-files", *boss_data_files,
           "--sampler", sampler,
           #"--sampler-config", "maxlike_tolerance", "1.0e-1",
           "--sampler-config", "emcee_walker", "72",
           "--sampler-config", "emcee_covariance_file", "data/KiDS1000/mocks/3x2pt_sample_logA_parameter_covariance.txt",
           "--overwrite"]
    subprocess.run(["python", script] + cmd, check=True)

    sampler = "multinest"
    run_name_root = "multinest_tol_0.01"
    run_name = f"{run_name_root}_{run_type}"
    cmd = ["--root-dir", root_dir,
           "--run-name", run_name,
           "--run-type", run_type,
           "--KiDS-data-file", twopoint_file,
           "--dz-covariance-file", dz_cov_file,
           "--BOSS-data-files", *boss_data_files,
           "--sampler", sampler,
           "--sampler-config", "nested_sampling_tolerance", "1.0e-2",
           "--overwrite"]
    subprocess.run(["python", script] + cmd, check=True)

    sampler = "multinest"
    run_name_root = "multinest_eff_0.1"
    run_name = f"{run_name_root}_{run_type}"
    cmd = ["--root-dir", root_dir,
           "--run-name", run_name,
           "--run-type", run_type,
           "--KiDS-data-file", twopoint_file,
           "--dz-covariance-file", dz_cov_file,
           "--BOSS-data-files", *boss_data_files,
           "--sampler", sampler,
           "--sampler-config", "multinest_efficiency", "0.1",
           "--overwrite"]
    subprocess.run(["python", script] + cmd, check=True)