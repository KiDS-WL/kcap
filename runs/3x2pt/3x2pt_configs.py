import subprocess

if __name__ == "__main__":
    script = "utils/run_kcap.py"
    root_dir = "runs/3x2pt/data/"
    
    # Cosmic shear + GGL twopoint file
    twopoint_file_template = "../Cat_to_Obs_K1000_P1/data/kids/fits/bp_KIDS1000_Blind{blind}_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits"

    # Covariance of the n(z)
    dz_cov_file = "../Cat_to_Obs_K1000_P1/data/kids/nofz/SOM_cov_multiplied.asc"

    # BOSS files
    boss_data_files = ["../Cat_to_Obs_K1000_P1/data/boss/Sanchez_etal_2017/BOSS.DR12.lowz.3xiwedges_measurements.txt",
                       "../Cat_to_Obs_K1000_P1/data/boss/Sanchez_etal_2017/BOSS.DR12.highz.3xiwedges_measurements.txt"]
    boss_cov_files =  ["../Cat_to_Obs_K1000_P1/data/boss/Sanchez_etal_2017/BOSS.DR12.lowz.3xiwedges_covmat.txt",
                       "../Cat_to_Obs_K1000_P1/data/boss/Sanchez_etal_2017/BOSS.DR12.highz.3xiwedges_covmat.txt"]


    multinest_settings = ["--sampler-config", "multinest_efficiency", "0.3",
                          "--sampler-config", "nested_sampling_tolerance", "1.0e-2",
                          "--sampler-config", "live_points", "250", # For final setup we probably want something higher than this
                         ]

    blinds = ["A",]                       # For final setup: ["A", "B", "C"]
    run_types = ["EE", "EE_w", "EE_nE_w"] # For final setup: ["EE", "nE", "w", "EE_nE", "EE_w", "nE_w", "EE_nE_w"]

    for blind in blinds:
        print(f"Blind {blind}")
        twopoint_file = twopoint_file_template.format(blind=blind)

        for run_type in run_types:
            print(f"  Run type: {run_type}")

            for sampler in ["test", "multinest"]:
                run_name_root = sampler

                run_name = f"{run_name_root}_blind{blind}_{run_type}"

                cmd = ["--root-dir", root_dir,
                        "--run-name", run_name,
                        "--run-type", run_type,
                        "--KiDS-data-file", twopoint_file,
                        "--dz-covariance-file", dz_cov_file,
                        "--BOSS-data-files", *boss_data_files,
                        "--BOSS-covariance-files", *boss_cov_files,
                        "--sampler", sampler,]

                if sampler == "multinest":
                    cmd += multinest_settings

                cmd += ["--overwrite"]

                subprocess.run(["python", script] + cmd, check=True)