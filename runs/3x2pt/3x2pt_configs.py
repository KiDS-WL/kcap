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
                          "--sampler-config", "live_points", "250", 
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


    # Cosmology chains
    root_dir = "runs/3x2pt/data_iterated_cov/cosmology/"

    blinds = []      
    run_types = ["EE", "EE_nE", "EE_w"]

    use_Planck = [False]

    for blind in blinds:
        print(f"Blind {blind}")
        twopoint_file = twopoint_file_template.format(blind=blind)

        for run_type in run_types:
            print(f"  Run type: {run_type}")

            for with_Planck in use_Planck:
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

                    if "w" in run_type:
                        # Allow wedges to time out
                        cmd += timeout_setttings

                    # cmd += ["--overwrite"]

                    subprocess.run(["python", script] + cmd, check=True)


    # Systematics chains
    root_dir = "runs/3x2pt/data_iterated_cov_very_fast/systematics/"

    # Configs for tomographic bin cuts
    tomographic_bin_cut_configs = []
    for cut_bin in [[0], [1], [2], [3], [4], [0,1]]:
        cut_tomographic_EE = [f"{i+1}+{j+1}" for i in range(5) for j in range(i,5) if i in cut_bin or j in cut_bin]
        cut_tomographic_nE = [f"{i+1}+{j+1}" for i in range(2) for j in range(5) if j in cut_bin]

        for extra_GGL_cut in ["1+1", "2+1", "2+2", "2+3"]:
            if extra_GGL_cut not in cut_tomographic_nE:
                cut_tomographic_nE.append(extra_GGL_cut)

        # Covariance of the n(z)
        cut_bin_str = '_'.join([str(c+1) for c in cut_bin])
        dz_cut_bin_cov_file = f"../Cat_to_Obs_K1000_P1/data/kids/nofz/SOM_cov_multiplied_bin{cut_bin_str}_removed.asc"
        
        # Compute the decorrelated dz mean shifts
        dz_cut_bin_cov = np.loadtxt(dz_cut_bin_cov_file)
        L_cut_bin = np.linalg.cholesky(dz_cut_bin_cov) 
        L_cut_bin_inv = np.linalg.inv(L_cut_bin)
        dx_cut_bin_mean = L_cut_bin_inv @ dz_mean

        parameter_settings = []
        for i, m in enumerate(dx_cut_bin_mean):
            if i in cut_bin:
                parameter_settings += ["--set-parameters", "nofz_shifts", f"p_{i+1}", f"{m}"]
            else:
                parameter_settings += ["--set-parameters", "nofz_shifts", f"p_{i+1}", f"-5.0 {m} 5.0"]
                parameter_settings += ["--set-priors", "nofz_shifts", f"p_{i+1}", f"gaussian {m} 1.0"]

        config = ["--set-keys", "scale_cuts", "cut_pair_PeeE", " ".join(cut_tomographic_EE),
                  "--set-keys", "scale_cuts", "cut_pair_PneE", " ".join(cut_tomographic_nE)]
        config += ["--dz-covariance-file", dz_cut_bin_cov_file]

        config += parameter_settings

        tomographic_bin_cut_configs.append((f"cut_z_bin_{''.join([str(c+1) for c in cut_bin])}", config))

    no_baryon_configs = [("no_baryon",   ["--set-parameters", "halo_model_parameters", "A", "3.13"]),]

    no_higher_order_configs = [# ("fix_ho_bias", ["--set-parameters", "bias_parameters", "b2_bin_1", "0.2",
                               #                  "--set-parameters", "bias_parameters", "gamma3_bin_1", "0.9",
                               #                  "--set-parameters", "bias_parameters", "a_vir_bin_1", "3.8",
                               #                  "--set-parameters", "bias_parameters", "b2_bin_2", "0.5",
                               #                  "--set-parameters", "bias_parameters", "gamma3_bin_2", "0.1",
                               #                  "--set-parameters", "bias_parameters", "a_vir_bin_2", "3.0",]),
                               ("zero_ho",      ["--set-keys", "wedges", "local_lag_g2", "F",
                                                 "--set-parameters", "bias_parameters", "b2_bin_1", "0.0",
                                                 "--set-parameters", "bias_parameters", "gamma2_bin_1", "0.0",
                                                 "--set-parameters", "bias_parameters", "gamma3_bin_1", "0.0",
                                                 "--set-parameters", "bias_parameters", "a_vir_bin_1", "0.0",
                                                 "--set-parameters", "bias_parameters", "b2_bin_2", "0.0",
                                                 "--set-parameters", "bias_parameters", "gamma3_bin_2", "0.0",
                                                 "--set-parameters", "bias_parameters", "gamma2_bin_2", "0.0",
                                                 "--set-parameters", "bias_parameters", "a_vir_bin_2", "0.0",])]
    
    configs = no_higher_order_configs + tomographic_bin_cut_configs + no_baryon_configs
    
    blinds = ["C"]
    for blind in blinds:
        print(f"Blind {blind}")
        twopoint_file = twopoint_file_template.format(blind=blind)

        run_type = "EE_nE_w"
        print(f"  Run type: {run_type}")

        for config_name, config in configs:
            print(f"    Config: {config_name}")
            for sampler in ["test", "multinest"]:
                run_name_root = sampler

                run_name = f"{run_name_root}_blind{blind}_{run_type}_{config_name}"

                # Base setup
                cmd = ["--root-dir", root_dir,
                        "--run-name", run_name,
                        "--run-type", run_type,
                        "--KiDS-data-file", twopoint_file,
                        "--BOSS-data-files", *boss_data_files,
                        "--BOSS-covariance-files", *boss_cov_files,
                        "--sampler", sampler,]      

                # GGL scale cuts
                if "nE" in run_type:
                    cmd += nE_scale_cuts

                if not "cut_z_bin" in config_name:
                    cmd += ["--dz-covariance-file", dz_cov_file]
                    # dz prior means
                    for i, m in enumerate(dx_mean):
                        cmd += ["--set-parameters", "nofz_shifts", f"p_{i+1}", f"-5.0 {m} 5.0"]
                        cmd += ["--set-priors", "nofz_shifts", f"p_{i+1}", f"gaussian {m} 1.0"]

                cmd += config

                # sampler settings
                if sampler == "multinest":
                    cmd += multinest_settings

                # Allow wedges to time out
                if "w" in run_type:
                    cmd += timeout_setttings

                # cmd += ["--overwrite"]

            subprocess.run(["python", script] + cmd, check=True)