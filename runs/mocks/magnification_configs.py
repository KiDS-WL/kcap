import subprocess

if __name__ == "__main__":

    # python utils/run_kcap.py --create-mocks --noiseless-mocks --root-dir data/KiDS1000/mocks/noisefree_lens_dz_0.01/ --KiDS-data-file ../K1000-data/Phase-1/twopoint/twoPoint_PneE+PeeE_mean_None_cov_theoryBuceros_nOfZ_bucerosBroad.fits --run-type EE_nE_w --run-name EE_nE_w_base --overwrite

    script = "utils/run_kcap.py"
    root_dir = "runs/mocks/magnification/"
    
    base_twopoint_file = "runs/mocks/data/noisefree_fiducials/base_EE_nE_w/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryBuceros_nOfZ_bucerosBroad_mock_noiseless.fits"
    
    # Need to change
    base_dz_cov_file = "data/KV450/nofz/DIR_cov.asc"

    base_boss_data_files = ["runs/mocks/data/noisefree_fiducials/base_EE_nE_w/data/BOSS/BOSS_mock_noiseless_bin_1.txt",
                            "runs/mocks/data/noisefree_fiducials/base_EE_nE_w/data/BOSS/BOSS_mock_noiseless_bin_2.txt"]
    
    # base_boss_cov_files = ["data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt",
    #                        "data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt"]

    # Main chains

    twopoint_file = base_twopoint_file
    dz_cov_file = base_dz_cov_file
    boss_data_files = base_boss_data_files

    sampler = "multinest"
    
    run_type = "EE_nE_w_magnification"
    run_name_root = "fast"
    for alpha in [1.0, 
                  3.0,
                  ]:
        run_name = f"{run_name_root}_{run_type}_alpha_{alpha}"
        cmd = ["--root-dir", root_dir,
               "--run-name", run_name,
               "--run-type", run_type,
               "--KiDS-data-file", twopoint_file,
               "--dz-covariance-file", dz_cov_file,
               "--BOSS-data-files", *boss_data_files,
               "--sampler", sampler,
               "--sampler-config", "multinest_efficiency", "0.3",
               "--sampler-config", "nested_sampling_tolerance", "1.0e-2",
               "--magnification-alpha", str(alpha), 
               "--overwrite"]
        subprocess.run(["python", script] + cmd, check=True)

    run_type = "EE_nE_w_magnification"
    run_name_root = "test_sampler"
    for alpha in [1.0, 
                  3.0,
                  ]:
        run_name = f"{run_name_root}_{run_type}_alpha_{alpha}"
        cmd = ["--root-dir", root_dir,
               "--run-name", run_name,
               "--run-type", run_type,
               "--KiDS-data-file", twopoint_file,
               "--dz-covariance-file", dz_cov_file,
               "--BOSS-data-files", *boss_data_files,
               "--sampler", "test",
               "--magnification-alpha", str(alpha), 
               "--overwrite"]
        subprocess.run(["python", script] + cmd, check=True)

