import subprocess

if __name__ == "__main__":

    # utils/run_kcap.py --create-mocks --noiseless-mocks --root-dir data/KiDS1000/mocks/noisefree/ --KiDS-data-file data/KiDS1000/mocks/twoPoint_PneE+PeeE_mean_theoryNoiseFree_cov_theoryBuceros_nOfZ_bucerosBroad.fits --run-type EE_nE_w --run-name EE_nE_w_base --overwrite

    script = "utils/run_kcap.py"
    sampler = "multinest"
    root_dir = "runs/methodology/chains"
    
    base_twopoint_file = "data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/KiDS/twoPoint_PneE+PeeE_mean_theoryNoiseFree_cov_theoryBuceros_nOfZ_bucerosBroad_mock_noiseless.fits"
    
    # Need to change
    base_dz_cov_file = "data/KV450/nofz/DIR_cov.asc"

    base_boss_data_files = ["data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS_mock_noiseless_bin_1.txt",
                            "data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS_mock_noiseless_bin_2.txt"]
    
    # base_boss_cov_files = ["data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt",
    #                        "data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt"]

    # Main chains
    run_name_root = "base"
    twopoint_file = base_twopoint_file
    dz_cov_file = base_dz_cov_file
    boss_data_files = base_boss_data_files
    for run_type in ["EE_nE_w",
                     "EE_nE",
                     "EE_w",
                     "nE_w",
                     "EE",
                     "nE",
                     "w"]:
        run_name = f"{run_name_root}_{run_type}"
        cmd = ["--root-dir", root_dir,
               "--run-name", run_name,
               "--run-type", run_type,
               "--KiDS-data-file", twopoint_file,
               "--dz-covariance-file", dz_cov_file,
               "--BOSS-data-files", *boss_data_files,
               "--sampler", sampler]
        subprocess.run(["python", script] + cmd, check=True)

    # Covariance chains
    run_name_root = "covariance"
    run_type = "EE_nE"
    dz_cov_file = base_dz_cov_file
    boss_data_files = base_boss_data_files
    # Loop over tuples of (name, twopoint-filename).
    for cov_name, twopoint_file in [("simple_theory", "data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/KiDS/twoPoint_PneE+PeeE_mean_theoryNoiseFree_cov_theoryBuceros_nOfZ_bucerosBroad_mock_noiseless.fits")]:
        run_name = f"{run_name_root}_{run_type}_{cov_name}"
        cmd = ["--root-dir", root_dir,
               "--run-name", run_name,
               "--run-type", run_type,
               "--KiDS-data-file", twopoint_file,
               "--dz-covariance-file", dz_cov_file,
               "--BOSS-data-files", *boss_data_files,
               "--fix-values", "nofz_shifts",     # Fix nuisance parameters.
               "--sampler", sampler]
        subprocess.run(["python", script] + cmd, check=True)

    # Prior chains
    run_name_root = "priors"
    run_type = "EE_nE"

    # Uncorrelated n(z) shifts
    twopoint_file = base_twopoint_file
    boss_data_files = base_boss_data_files
    # File with uncorrorelated dz prior covariance
    dz_cov_file = "data/KV450/nofz/id_cov.asc"
    prior_name = "uncorrelated_nz_shifts"
    run_name = f"{run_name_root}_{run_type}_{prior_name}"
    cmd = ["--root-dir", root_dir,
            "--run-name", run_name,
            "--run-type", run_type,
            "--KiDS-data-file", twopoint_file,
            "--dz-covariance-file", dz_cov_file,
            "--BOSS-data-files", *boss_data_files,
            "--sampler", sampler]
    subprocess.run(["python", script] + cmd, check=True)

    # bin, scale cut chains
    run_name_root = "cuts"
    twopoint_file = base_twopoint_file
    dz_cov_file = base_dz_cov_file
    boss_data_files = base_boss_data_files
    for cut_name, set_keys_args in [("cut_overlap_ggl", ("scale_cuts", "cut_pair_PneE", "1+1 1+2 1+3 2+1 2+2 2+3 2+4")),
                                    ("conservative_ell_cut_ggl", ("scale_cuts", "keep_ang_PneE", "100 500"))]:
        run_name = f"{run_name_root}_{run_type}_{cut_name}"
        cmd = ["--root-dir", root_dir,
               "--run-name", run_name,
               "--run-type", run_type,
               "--KiDS-data-file", twopoint_file,
               "--dz-covariance-file", dz_cov_file,
               "--BOSS-data-files", *boss_data_files,
               "--set-keys", *set_keys_args,  
               "--sampler", sampler]
        subprocess.run(["python", script] + cmd, check=True)