import subprocess
import numpy as np

if __name__ == "__main__":
    script = "utils/run_kcap.py"
    
    # Cosmic shear + GGL twopoint file
    twopoint_file_template = "../Cat_to_Obs_K1000_P1/data/kids/fits_iterative_covariance/bp_KIDS1000_Blind{blind}_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits"

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
                          "--sampler-config", "live_points", "500", 
                         ]

    timeout_setttings = ["--set-keys", "wedges", "timeout", "600.0"]

    nE_scale_cuts = ["--set-keys", "scale_cuts", "keep_ang_PneE_1_1", "100 300",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_1_2", "100 300",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_1_3", "100 300",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_1_4", "100 300",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_1_5", "100 300",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_2_1", "100 600",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_2_2", "100 600",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_2_3", "100 600",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_2_4", "100 600",
                     "--set-keys", "scale_cuts", "keep_ang_PneE_2_5", "100 600",]

    Planck_settings = ["--enable-modules", "planck_like",
                       "--set-keys", "camb", "mode", "cmb",
                       "--set-keys", "camb", "lmax", "2650",
                       "--set-keys", "camb", "nonlinear", "both",
                       "--set-keys", "camb", "do_lensing", "T",
                       "--set-keys", "camb", "do_reionization", "T",
                       # Planck TTTEEE+lowl+lowE 5 sigma ranges, with S8 and ns having a 7 sigma lower range and h having a 7 sigma upper range.
                       "--set-parameters", "cosmological_parameters", "omch2",     "0.11336956221293837   0.12       0.12703091064106695",
                       "--set-parameters", "cosmological_parameters", "ombh2",     "0.021615770756319073  0.0225     0.023103729722525095",
                       "--set-parameters", "cosmological_parameters", "h0",        "0.6425290065540109    0.7        0.7150188763839126",
                       "--set-parameters", "cosmological_parameters", "n_s",        "0.9342806461291905    0.97       0.986698010466018",
                       "--set-parameters", "cosmological_parameters", "S_8_input", "0.7229290272333152    0.7458     0.9134139168266714",
                       "--set-parameters", "cosmological_parameters", "tau",       "0.015070795999054837  0.0543     0.09381939280201967",
                       "--set-parameters", "planck",                  "a_planck",  "0.9879083109867925    1.000610   1.0130810744845216",
                       "--set-priors", "planck", "a_planck", "gaussian 1.0 0.0025",]

    Planck_lensing_settings = ["--enable-modules", "planck_lensing_like",
                               "--set-keys", "camb", "mode", "cmb",
                               "--set-keys", "camb", "lmax", "2650",
                               "--set-keys", "camb", "nonlinear", "both",
                               "--set-keys", "camb", "do_lensing", "T",
                               "--set-keys", "camb", "do_reionization", "T",
                               "--set-keys", "camb", "lens_potential_accuracy", "1",
                               "--set-parameters", "cosmological_parameters", "tau",       "0.055",]

    Pantheon_settings = ["--enable-modules", "pantheon_like",
                         "--set-parameters", "supernova_params",  "M",  "-22.0  -19.0   -17.0",]

    Pantheon_fixed_cal_settings = ["--enable-modules", "pantheon_like",
                                   "--set-parameters", "supernova_params",  "M",  "-19.0",]

    Riess2019_settings = ["--set-priors", "cosmological_parameters", "h0",  "gaussian 0.7403 0.0142"]

    Pantheon_configs = [("LCDM_Pantheon",  Pantheon_settings),]
    Pantheon_fixed_cal_configs = [("LCDM_Pantheon_fixed_cal",  Pantheon_fixed_cal_settings),]
    Riess2019_configs = [("LCDM_Riess2019",  Riess2019_settings),]

    Planck_lensing_configs = [("LCDM_Planck_lensing", Planck_lensing_settings)]

    # Cosmology chains
    root_dir = "runs/extended_cosmologies/data/cosmology/"

    blinds = ["C"]      
    run_types = ["EE_nE_w"]

    use_Planck = [False]

    configs = Riess2019_configs

    for blind in blinds:
        print(f"Blind {blind}")
        twopoint_file = twopoint_file_template.format(blind=blind)

        for run_type in run_types:
            print(f"  Run type: {run_type}")

            for with_Planck in use_Planck:
                print(f"    Include Planck: {with_Planck}")

                for config_name, config in configs:
                    print(f"      Config: {config_name}")
                    for sampler in ["test", "multinest"]:
                        run_name_root = sampler

                        if sampler == "multinest":
                            if multinest_settings[-1] == "500":
                                run_name_root += "_medium"
                            elif multinest_settings[-1] == "250":
                                run_name_root += "_fast"

                        run_name = f"{run_name_root}_blind{blind}_{run_type}_{config_name}"
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

                        cmd += config
                        
                        # sampler settings
                        if sampler == "multinest":
                            cmd += multinest_settings

                        if "w" in run_type:
                            # Allow wedges to time out
                            cmd += timeout_setttings

                        # cmd += ["--overwrite"]

                        subprocess.run(["python", script] + cmd, check=True)
