import subprocess

if __name__ == "__main__":
    script = "utils/run_kcap.py"
    
    base_twopoint_file = "runs/extended_cosmologies/mocks/data/noisefree/lcdm_EE_nE_w/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noiseless.fits"
    
    base_dz_cov_file = "data/KV450/nofz/SOM_cov.asc"

    base_boss_data_files = ["runs/extended_cosmologies/mocks/data/noisefree/lcdm_EE_nE_w/data/BOSS/BOSS_mock_noiseless_bin_1.txt",
                            "runs/extended_cosmologies/mocks/data/noisefree/lcdm_EE_nE_w/data/BOSS/BOSS_mock_noiseless_bin_2.txt"]
    
    # base_boss_cov_files = ["data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS.DR12.lowz.3xiwedges_covmat.txt",
    #                        "data/KiDS1000/mocks/noisefree/EE_nE_w_base/data/BOSS/BOSS.DR12.highz.3xiwedges_covmat.txt"]

    # Main chains

    root_dir = "runs/extended_cosmologies/mocks/lcdm_mock_fast_MG_sigma8/"

    twopoint_file = base_twopoint_file
    dz_cov_file = base_dz_cov_file
    boss_data_files = base_boss_data_files

    for sampler in ["test", "multinest"]:

        run_type = "EE_fR"
        run_name_root = sampler

        for name, extra_parameter in {#"LCDM"           : [("cosmological_parameters", "log10_fR0", "-12.0")],
                                      #"fRCDM"          : [("cosmological_parameters", "fR0", "1e-8 1e-7 1e-4")],
                                      "fRCDM_logfR0_8_2"   : ["--set-parameters", "cosmological_parameters", "log10_fR0", "-8.0 -7.0 -2.0"],
                                      #"fRCDM_logfR0_12_3"   : [("cosmological_parameters", "log10_fR0", "-12.0 -7.0 -3.0"),
                                      #"fRCDM_logfR0_12_3_wrap"   : [("cosmological_parameters", "log10_fR0", "-12.0 -7.0 -3.0")],
                                      #"fRCDM_logfR0_8_3_wrap"   : [("cosmological_parameters", "log10_fR0", "-8.0 -7.0 -3.0")],
                                    #   "fRCDM_logfR0_8_3_folded"   : ["--set-parameters", "cosmological_parameters", "log10_fR0_folded", "-13.0 -8.0 -3.0",
                                    #                                  "--enable-module", "sample_folded_prior",
                                    #                                  "--set-keys", "sample_folded_prior", "name", "cosmological_parameters/log10_fR0",
                                    #                                  "--set-keys", "sample_folded_prior", "fold", "-8.0",],
                                      # "fRCDM_logfR0_8_2_folded"   : ["--set-parameters", "cosmological_parameters", "log10_fR0_folded", "-14.0 -8.0 -2.0",
                                      #                                "--enable-module", "sample_folded_prior",
                                      #                                "--set-keys", "sample_folded_prior", "name", "cosmological_parameters/log10_fR0",
                                      #                                "--set-keys", "sample_folded_prior", "fold", "-8.0",],
                                      "fRCDM_logfR0_8_2_wrapped"   : ["--set-parameters", "cosmological_parameters", "log10_fR0", "-8.0 -7.0 -2.0",
                                                                      "--sampler-config", "wrapped_params", "cosmological_parameters--log10_fr0",],
                                      #"fRCDM_logfR0_8_2_emcee"   : ["--set-parameters", "cosmological_parameters", "log10_fR0", "-8.0 -5.0 -2.0",],
                                    #   "fRCDM_fR0_0_0p001_folded"  : ["--set-parameters", "cosmological_parameters", "fR0_folded", "-1e-3 0.0 1e-3",
                                    #                                  "--enable-module", "sample_folded_prior",
                                    #                                  "--set-keys", "sample_folded_prior", "name", "cosmological_parameters/fR0",
                                    #                                  "--set-keys", "sample_folded_prior", "fold", "0.0",
                                    #                                  "--set-keys", "reaction", "log10_fR0", "F",],
                                        }.items():
            run_name = f"{run_name_root}_{name}_{run_type}"
            cmd = ["--root-dir", root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--BOSS-data-files", *boss_data_files,
                    "--sampler", sampler,]
            
            cmd += extra_parameter

            # Adjust priors
            cmd += ["--set-parameters", "cosmological_parameters", "S_8_input", "0.45 0.7458 1.3"]

            cmd += ["--derived-parameters", "cosmological_parameters/sigma_8_LCDM", "cosmological_parameters/S_8_LCDM"]

            # Speed things up
            # cmd += ["--set-keys", "camb", "z_mid", "2.0"]
            # cmd += ["--set-keys", "camb", "nz_mid", "50"]
            # cmd += ["--set-keys", "camb", "nz", "100"]
            
            cmd += ["--set-keys", "reaction", "mode", "f(R)"]

            if sampler == "multinest":
                cmd += ["--sampler-config", "multinest_efficiency", "0.3",
                        "--sampler-config", "nested_sampling_tolerance", "1.0e-2",
                        "--sampler-config", "live_points", "500",
                        ]

            cmd += ["--overwrite"]

            subprocess.run(["python", script] + cmd, check=True)