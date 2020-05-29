import subprocess
import numpy as np

if __name__ == "__main__":
    script = "utils/run_kcap.py"
    root_dir = "runs/3x2pt/data_w_dz_prior_means/"
    
    # Cosmic shear + GGL twopoint file
    twopoint_file_template = "../Cat_to_Obs_K1000_P1/data/kids/fits/bp_KIDS1000_Blind{blind}_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits"

    # Covariance of the n(z)
    dz_cov_file = "../Cat_to_Obs_K1000_P1/data/kids/nofz/SOM_cov_multiplied.asc"
    dz_mean_file = "../Cat_to_Obs_K1000_P1/data/kids/nofz/deltaz.asc"
    
    # Compute the decorrelated dz mean shifts
    dz_cov = np.loadtxt(dz_cov_file)
    dz_mean = np.loadtxt(dz_mean_file)
    L = np.linalg.cholesky(dz_cov) 
    L_inv = np.linalg.inv(L)
    dx_mean = L_inv @ dz_mean

    # BOSS files
    boss_data_files = ["../Cat_to_Obs_K1000_P1/data/boss/Sanchez_etal_2017/BOSS.DR12.lowz.3xiwedges_measurements.txt",
                       "../Cat_to_Obs_K1000_P1/data/boss/Sanchez_etal_2017/BOSS.DR12.highz.3xiwedges_measurements.txt"]
    boss_cov_files =  ["../Cat_to_Obs_K1000_P1/data/boss/Sanchez_etal_2017/BOSS.DR12.lowz.3xiwedges_covmat.txt",
                       "../Cat_to_Obs_K1000_P1/data/boss/Sanchez_etal_2017/BOSS.DR12.highz.3xiwedges_covmat.txt"]


    multinest_settings = ["--sampler-config", "multinest_efficiency", "0.3",
                          "--sampler-config", "nested_sampling_tolerance", "1.0e-2",
                          "--sampler-config", "live_points", "250", # For final setup we probably want something higher than this
                         ]

    nE_scale_cuts = ["--set-keys", "scale_cuts", "keep_ang_PneE_1_2", "100 700",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_1_3", "100 700",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_1_4", "100 700",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_1_5", "100 700",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_2_4", "100 300",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_2_5", "100 300",]

    Planck_settings = ["--enable-modules", "planck_like",
                       "--set-keys", "camb", "mode", "cmb",
                       "--set-keys", "camb", "lmax", "2650",
                       "--set-keys", "camb", "nonlinear", "both",
                       "--set-keys", "camb", "do_lensing", "T",
                       "--set-keys", "camb", "do_reionization", "T",
                       "--set-parameters", "cosmological_parameters", "tau", "0.015070795999054837    0.0543     0.09381939280201967",
                       "--set-parameters", "planck", "a_planck", "0.9879083109867925 1.000610 1.0130810744845216",
                       "--set-priors", "planck", "a_planck", "gaussian 1.0 0.0025",]


    blinds = ["A",]                       # For final setup: ["A", "B", "C"]
    run_types = ["EE", "EE_w", "EE_nE_w", "w"] # For final setup: ["EE", "nE", "w", "EE_nE", "EE_w", "nE_w", "EE_nE_w"]

    for blind in blinds:
        print(f"Blind {blind}")
        twopoint_file = twopoint_file_template.format(blind=blind)

        for run_type in run_types:
            print(f"  Run type: {run_type}")

            for with_Planck in [False, True]:
                print(f"    Include Planck: {with_Planck}")

                for sampler in ["test", "multinest"]:
                    run_name_root = sampler

                    run_name = f"{run_name_root}_blind{blind}_{run_type}"
                    if with_Planck:
                        run_name += "_Planck"

                    # Base setup
                    cmd = ["--root-dir", root_dir,
                            "--run-name", run_name,
                            "--run-type", run_type,
                            "--KiDS-data-file", twopoint_file,
                            "--dz-covariance-file", dz_cov_file,
                            "--BOSS-data-files", *boss_data_files,
                            "--BOSS-covariance-files", *boss_cov_files,
                            "--sampler", sampler,]

                    # dz prior means
                    for i, m in enumerate(dx_mean):
                        cmd += ["--set-parameters", "nofz_shifts", f"p_{i+1}", f"-5.0 {m} 5.0"]
                        cmd += ["--set-priors", "nofz_shifts", f"p_{i+1}", f"gaussian {m} 1.0"]

                    # GGL scale cuts
                    if "nE" in run_type:
                        cmd += nE_scale_cuts

                    # Add Planck likelihood
                    if with_Planck:
                        cmd += Planck_settings

                    # sampler settings
                    if sampler == "multinest":
                        cmd += multinest_settings

                    cmd += ["--overwrite"]

                    subprocess.run(["python", script] + cmd, check=True)