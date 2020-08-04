import subprocess
import numpy as np

if __name__ == "__main__":
    script = "utils/run_kcap.py"
    
    # Cosmic shear + GGL twopoint file
    twopoint_file_template = "../Cat_to_Obs_K1000_P1/data/kids/fits/bp_KIDS1000_Blind{blind}_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits"

    # Covariance of the n(z)
    dz_cov_file = "../Cat_to_Obs_K1000_P1/data/kids/nofz/SOM_cov_multiplied.asc"
    dz_mean_file = "../Cat_to_Obs_K1000_P1/data/kids/nofz/deltaz.asc"

    # BOSS files
    boss_data_files = ["../Cat_to_Obs_K1000_P1/data/boss/Sanchez_etal_2017/BOSS.DR12.lowz.3xiwedges_measurements.txt",
                       "../Cat_to_Obs_K1000_P1/data/boss/Sanchez_etal_2017/BOSS.DR12.highz.3xiwedges_measurements.txt"]
    boss_cov_files =  ["../Cat_to_Obs_K1000_P1/data/boss/Sanchez_etal_2017/BOSS.DR12.lowz.3xiwedges_covmat.txt",
                       "../Cat_to_Obs_K1000_P1/data/boss/Sanchez_etal_2017/BOSS.DR12.highz.3xiwedges_covmat.txt"]


    multinest_settings = ["--sampler-config", "multinest_efficiency", "0.3",
                          "--sampler-config", "nested_sampling_tolerance", "1.0e-2",
                          "--sampler-config", "live_points", "250", # For final setup we probably want something higher than this
                         ]

    timeout_setttings = ["--set-keys", "wedges", "timeout", "600.0"]
    only_linear_pofk_settings = ["--set-keys", "camb", "nonlinear", "none"]

    # Cosmology chains
    root_dir = "runs/3x2pt/clustering_only_fast/cosmology/"

    run_type = "w"

    twopoint_file = twopoint_file_template.format(blind="A")

    for sampling_choice in ["fiducial", "S8_wide_ns"]:
        for sampler in ["test", "multinest"]:
            run_name_root = sampler

            run_name = f"{run_name_root}_sample_{sampling_choice}_{run_type}"

            # Base setup
            cmd = ["--root-dir", root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--BOSS-data-files", *boss_data_files,
                    "--BOSS-covariance-files", *boss_cov_files,
                    "--sampler", sampler,]

            # sampler settings
            if sampler == "multinest":
                cmd += multinest_settings

            # Allow wedges to time out
            cmd += timeout_setttings
            # Don't need HMCode non-linear P(k)
            cmd += only_linear_pofk_settings

            if sampling_choice == "S8_wide_ns":
                cmd += ["--set-parameters", "cosmological_parameters", "n_s", "0.5 0.97 1.1"]

            if sampling_choice == "logAs":
                cmd += ["--enable-modules", "sample_ln_As"]
                cmd += ["--cut-modules", "sample_S8",]
                cmd += ["--cut-modules", "sigma8toAs",]
                cmd += ["--set-parameters", "cosmological_parameters", "ln_1e10_a_s", "1.5 2.72 4.0"]
                cmd += ["--set-parameters", "cosmological_parameters", "S_8_input", "0.7458"]


            #cmd += ["--overwrite"]

            subprocess.run(["python", script] + cmd, check=True)