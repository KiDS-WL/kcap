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

        for run_type in ["EE"]:#, "EE_nE_w", "w"]:
            run_name_root = "fast" if sampler == "multinest" else "test"

            for name, extra_parameter in {#"LCDM_re"   : [],
                                          # "kCDM"   : [("cosmological_parameters", "omega_k", "-0.3 0.0 0.3")],
                                           "nuCDM_sym_prior"  : [("cosmological_parameters", "mnu_proxy", "-3.0 0.06 3.0")],
                                           #"wCDM"   : [("cosmological_parameters", "w", "-3.0 -1.0 -0.33")],
                                           #"waCDM"  : [("cosmological_parameters", "w", "-3.0 -1.0 -0.33"),
                                           #            ("cosmological_parameters", "wa", "-3.0 0.0 3.0")],
                                          }.items():
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
                    if p[1] == "mnu_proxy":
                        cmd += ["--set-parameters", "cosmological_parameters", "mnu", "none"]
                        cmd += ["--enable-modules", "sample_negative_mnu"]

                    if p[1] == "wa":
                        cmd += ["--set-keys", "camb", "use_ppf_w", "T"]

                if sampler == "multinest":
                    cmd += ["--sampler-config", "multinest_efficiency", "0.3",
                            "--sampler-config", "nested_sampling_tolerance", "1.0e-2",
                            #"--sampler-config", "live_points", "500",
                            #"--sampler-config", "multinest_wrapped_params", "cosmological_parameters--mnu",
                           ]

                cmd += ["--overwrite"]

                subprocess.run(["python", script] + cmd, check=True)