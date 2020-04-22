import subprocess

if __name__ == "__main__":
    script = "utils/run_kcap.py"
    root_dir = "runs/extended_cosmologies/mocks/"
    
    base_twopoint_file = "runs/extended_cosmologies/mocks/data/noisefree/lcdm_EE_nE_w/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits"
    
    base_dz_cov_file = "data/KV450/nofz/SOM_cov.asc"

    base_boss_data_files = ["runs/extended_cosmologies/mocks/data/noisefree/lcdm_EE_nE_w/data/BOSS/BOSS_mock_noiseless_bin_1.txt",
                            "runs/extended_cosmologies/mocks/data/noisefree/lcdm_EE_nE_w/data/BOSS/BOSS_mock_noiseless_bin_2.txt"]
    
    # base_boss_cov_files = ["data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt",
    #                        "data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt"]

    # Main chains

    twopoint_file = base_twopoint_file
    dz_cov_file = base_dz_cov_file
    boss_data_files = base_boss_data_files

    for sampler in ["test", "multinest"]:

        for run_type in ["EE", "EE_nE_w"]:
            run_name_root = "fast" if sampler == "multinest" else "test"

            for name, extra_parameter in {"kCDM"   : [("cosmological_parameters", "omega_k", "-0.2 0.0 0.2")],
                                          "nuCDM"  : [("cosmological_parameters", "mnu", "0.0 0.06 1.0")],
                                          "wCDM"   : [("cosmological_parameters", "w", "0.0 -1.0 -3.0")],
                                          "waCDM"  : [("cosmological_parameters", "w", "0.0 -1.0 -3.0"),
                                                      ("cosmological_parameters", "wa", "-0.5 0.0 0.5")],}.items():
                run_name = f"{run_name_root}_{name}_{run_type}"
                cmd = ["--root-dir", root_dir,
                        "--run-name", run_name,
                        "--run-type", run_type,
                        "--KiDS-data-file", twopoint_file,
                        "--dz-covariance-file", dz_cov_file,
                        "--BOSS-data-files", *boss_data_files,
                        "--sampler", sampler,]
                for p in extra_parameter:
                    cmd += ["--set-parameters", *p]

                if sampler == "multinest":
                    cmd += ["--sampler-config", "multinest_efficiency", "0.3",
                            "--sampler-config", "nested_sampling_tolerance", "1.0e-2",]

                cmd += ["--overwrite"]

                subprocess.run(["python", script] + cmd, check=True)