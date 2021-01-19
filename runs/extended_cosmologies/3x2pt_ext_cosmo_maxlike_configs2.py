
import subprocess
import os

import numpy as np

import sys
sys.path.append("utils")
import process_chains


if __name__ == "__main__":
    script = "utils/run_kcap.py"

    multinest_chain_files = {"wCDM" : {"EE_nE_w"        : "runs/extended_cosmologies/data/cosmology/multinest_blindC_EE_nE_w_wLCDM/chain/samples_multinest_blindC_EE_nE_w_wLCDM.txt",
                                       "EE"             : "runs/extended_cosmologies/data/cosmology/multinest_blindC_EE_wLCDM/chain/samples_multinest_blindC_EE_wLCDM.txt",
                                       "w"              : "runs/extended_cosmologies/data/cosmology/multinest_medium_blindC_w_wLCDM/chain/samples_multinest_medium_blindC_w_wLCDM.txt",
                                       "EE_nE_w_Planck" : "runs/extended_cosmologies/data/cosmology/multinest_medium_blindC_EE_nE_w_wCDM_Planck/chain/samples_multinest_medium_blindC_EE_nE_w_wCDM_Planck.txt",
                                       "Planck"         : "runs/3x2pt/Planck/multinest_Planck_wCDM/chain/samples_multinest_Planck_wCDM.txt",},
                             "oCDM" : {"EE_nE_w"        : "runs/extended_cosmologies/data/cosmology/multinest_blindC_EE_nE_w_oLCDM/chain/samples_multinest_blindC_EE_nE_w_oLCDM.txt",
                                       "EE"             : "runs/extended_cosmologies/data/cosmology/multinest_blindC_EE_oLCDM/chain/samples_multinest_blindC_EE_oLCDM.txt",
                                       "w"              : "runs/extended_cosmologies/data/cosmology/multinest_medium_blindC_w_oLCDM/chain/samples_multinest_medium_blindC_w_oLCDM.txt",
                                       "Planck"         : "runs/3x2pt/Planck/multinest_Planck_oCDM/chain/samples_multinest_Planck_oCDM.txt",
                                       },
                             "nuCDM": {"EE_nE_w"        : "runs/extended_cosmologies/data/cosmology/multinest_blindC_EE_nE_w_nuLCDM/chain/samples_multinest_blindC_EE_nE_w_nuLCDM.txt",
                                       "w"              : "runs/extended_cosmologies/data/cosmology/multinest_medium_blindC_w_nuLCDM/chain/samples_multinest_medium_blindC_w_nuLCDM.txt",
                                       "EE"             : "runs/extended_cosmologies/data/cosmology/multinest_blindC_EE_nuLCDM_wrap/chain/samples_multinest_blindC_EE_nuLCDM_wrap.txt",},
                             "fRCDM" :{ "EE_fR"         : "runs/extended_cosmologies/data/cosmology/multinest_medium_blindC_EE_fR_fRCDM/chain/samples_multinest_medium_blindC_EE_fR_fRCDM.txt",
                                       },
                             "fixAs" :{"EE"             : "runs/3x2pt/data_iterated_cov_fast_extra/systematics/multinest_blindC_EE_fix_As/chain/samples_multinest_blindC_EE_fix_As.txt",
                                       "EE_nE_w"        : "runs/3x2pt/data_iterated_cov_fast_extra/systematics/multinest_blindC_EE_nE_w_fix_As/chain/samples_multinest_blindC_EE_nE_w_fix_As.txt"},
                             "Pantheon" : {"EE"         : "runs/extended_cosmologies/data/cosmology/multinest_blindC_EE_LCDM_Pantheon/chain/samples_multinest_blindC_EE_LCDM_Pantheon.txt",
                                           "EE_nE_w"    : "runs/extended_cosmologies/data/cosmology/multinest_medium_blindC_EE_nE_w_LCDM_Pantheon/chain/samples_multinest_medium_blindC_EE_nE_w_LCDM_Pantheon.txt"},
                             "CMBlensing" : {"EE_nE_w"  : "runs/extended_cosmologies/data/cosmology/multinest_blindC_EE_nE_w_LCDM_Planck_lensing/chain/samples_multinest_blindC_EE_nE_w_LCDM_Planck_lensing.txt",
                                             "EE"       : "runs/extended_cosmologies/data/cosmology/multinest_blindC_EE_LCDM_Planck_lensing/chain/samples_multinest_blindC_EE_LCDM_Planck_lensing.txt"},
                             "HMCode2020" : {"EE"       : "runs/extended_cosmologies/data/cosmology/multinest_medium_blindC_EE_HMCode2020_wide/chain/samples_multinest_medium_blindC_EE_HMCode2020_wide.txt",
                                             "EE_nE_w"  : "runs/extended_cosmologies/data/cosmology/multinest_medium_blindC_EE_nE_w_HMCode2020_wide/chain/samples_multinest_medium_blindC_EE_nE_w_HMCode2020_wide.txt"},
                            }
    n_start_point = 36

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

    timeout_setttings = ["--set-keys", "wedges", "timeout", "600.0"]

    Planck_settings = ["--enable-modules", "planck_like",
                       "--set-keys", "camb", "mode", "cmb",
                       "--set-keys", "camb", "lmax", "2650",
                       "--set-keys", "camb", "nonlinear", "both",
                       "--set-keys", "camb", "do_lensing", "T",
                       "--set-keys", "camb", "do_reionization", "T",]

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

    wLCDM_configs = [("wCDM",          []),]
    oLCDM_configs = [("oCDM",          []),]
    nuLCDM_configs = [("nuCDM",          []),]
    fRCDM_configs = [("fRCDM",         ["--set-keys", "reaction", "mode", "f(R)",
                                        "--derived-parameters", "cosmological_parameters/sigma_8_LCDM", "cosmological_parameters/S_8_LCDM",]),]

    fixAs_configs = [("fixAs",          ["--set-parameters", "cosmological_parameters", "A_s", "2.1e-9",
                                         "--set-parameters", "cosmological_parameters", "S_8_input", "none",
                                         "--cut-modules", "sample_S8",
                                         "--cut-modules", "sigma8toAs"]),]

    CMBlensing_configs = [("CMBlensing",  ["--enable-modules", "planck_lensing_like",
                                           "--set-keys", "camb", "mode", "cmb",
                                           "--set-keys", "camb", "lmax", "2650",
                                           "--set-keys", "camb", "nonlinear", "both",
                                           "--set-keys", "camb", "do_lensing", "T",
                                           "--set-keys", "camb", "do_reionization", "T",
                                           "--set-keys", "camb", "lens_potential_accuracy", "1",
                                           "--set-parameters", "cosmological_parameters", "tau",       "0.055",]),]

    HMCode_configs = [("HMCode2020",  ["--set-parameters", "halo_model_parameters", "A", "none",
                                         "--cut-modules", "one_parameter_hmcode",
                                         "--set-keys", "camb", "halofit_version", "mead2020_feedback"])]

    Pantheon_configs = [("Pantheon",  ["--enable-modules", "pantheon_like",]),]

    configs = nuLCDM_configs

    blind = "C"
    run_types = [#"EE_nE_w", "EE",# "w", #"EE", "EE_nE",
                 #"EE_nE_w_Planck", 
                 #"Planck"
                 #"EE_fR"
                 "w"
                 ]

    sampler = "maxlike"
    run_name_root = "MAP"

    solver = "Nelder-Mead"

    for config_name, config in configs:
        print(f"Config {config_name}")
        root_dir = f"runs/extended_cosmologies/data_MAP/cosmology_{config_name}_NM/"

        twopoint_file = twopoint_file_template.format(blind=blind)
        multinest_chain_file = multinest_chain_files[config_name]

        for run_type in run_types:
            if run_type not in multinest_chain_file:
                continue

            print(f"  Run type: {run_type}")

            # Get starting points
            chain = process_chains.load_chain(multinest_chain_file[run_type], burn_in=0,
                                              keep_raw_samples=True, ignore_inf=True)

            param_names = [n.name for n in chain.getParamNames().names]
            #print(param_names)

            d = chain.getParams()
            max_logpost_idx = (d.logpost).argmax()
            logpost_sort_idx = np.argsort(d.logpost)[::-1]
            
            start_points = []
            for idx in logpost_sort_idx[:n_start_point]:
                start_points.append({n : getattr(d, n)[idx] for n in param_names})

                # Hacky due to getdist deleting parameters that are not changing
                if "EE" in run_type and "a_baryon" not in param_names and config_name != "HMCode2020" :
                    a_baryon = chain.raw_chain_samples[0,chain.raw_chain_parameter_names.index("a_baryon")]
                    start_points[-1]["a_baryon"] = a_baryon

            # if solver  == "Powell":
            #     with open(multinest_chain_file[run_type], "r") as f:
            #         for l in f.readlines():
            #             if "n_varied" in l:
            #                 n_varied = int(l.split("=")[-1])
            #                 break
            #     cov = chain.cov(params=list(range(n_varied)))
            #     v, e = np.linalg.eig(cov)

            #     step_size = 0.1
            #     search_directions = np.array([step_size * v[i] * e[:,i] for i in range(n_varied)])

            for i, p in enumerate(start_points):
                print(f"Starting point {i} logpost : {p['logpost']:.2f}", end=" ")
                try:
                    print(f"logprior: {p['logprior']:.2f}", end=" ")
                except:
                    pass
                print(f"loglike: {chain.loglikes[logpost_sort_idx[i]]:.2f}")

                run_name = f"{run_name_root}_{i}_blind{blind}_{run_type}_{solver}"

                # np.savetxt(os.path.join(root_dir, run_name, "config", "powell_search_directions.txt"), search_directions)

                # Base setup
                cmd = [ "--root-dir", root_dir,
                        "--run-name", run_name,
                        "--KiDS-data-file", twopoint_file,
                        "--dz-covariance-file", dz_cov_file,
                        "--BOSS-data-files", *boss_data_files,
                        "--BOSS-covariance-files", *boss_cov_files,
                        "--sampler", sampler,]

                if "Planck" in run_type:
                    cmd += Planck_settings
                    if config_name == "wCDM":
                        cmd += [# Planck TTTEEE+lowl+lowE wCDM 5 sigma ranges
                                "--set-parameters", "cosmological_parameters", "omch2",     f"0.11311880756197211  {p['omegach2']} 0.12673547660171416",
                                "--set-parameters", "cosmological_parameters", "ombh2",     f"0.021654218587807587 {p['omegabh2']} 0.023132100816368375",
                                "--set-parameters", "cosmological_parameters", "h0",        f"0.64                 {p['h']}        0.82",
                                "--set-parameters", "cosmological_parameters", "n_s",       f"0.9437568113540074   {p['ns']} 0.9870276396097222",
                                "--set-parameters", "cosmological_parameters", "S_8_input", f"0.6280620828064101   {p['s8proxy']} 0.9269044096801944",
                                "--set-parameters", "cosmological_parameters", "tau",       f"0.014360345677643029 {p['tau']} 0.09350636226730699",
                                "--set-parameters", "planck",                  "a_planck",  f"0.9881436963908558   {p['calPlanck']} 1.0129303406710994",
                                "--set-priors", "planck", "a_planck", "gaussian 1.0 0.0025",]
                    elif config_name == "oCDM":
                        cmd += [# Planck TTTEEE+lowl+lowE oCDM 5 sigma ranges, with S8 and ns having a 7 sigma lower range and omegak having a 7 sigma range.
                                "--set-parameters", "cosmological_parameters", "omch2",     f"0.11060872577816308  {p['omegach2']} 0.125577998484406",
                                "--set-parameters", "cosmological_parameters", "ombh2",     f"0.021743311058688833 {p['omegabh2']} 0.023454890484641437",
                                "--set-parameters", "cosmological_parameters", "h0",        f"0.64                 {p['h']}        0.82",
                                "--set-parameters", "cosmological_parameters", "n_s",       f"0.9369207025250993   {p['ns']} 0.9946939190599493",
                                "--set-parameters", "cosmological_parameters", "S_8_input", f"0.6408512866385077   {p['s8proxy']} 1.2236067309574206",
                                "--set-parameters", "cosmological_parameters", "tau",       f"0.007601071445446535 {p['tau']} 0.08967168037374854",
                                "--set-parameters", "planck",                  "a_planck",  f"0.9874241706445309   {p['calPlanck']} 1.0124671150845403",
                                "--set-priors", "planck", "a_planck", "gaussian 1.0 0.0025",]

                    else:
                        cmd += [# Planck TTTEEE+lowl+lowE 5 sigma ranges, with S8 and ns having a 7 sigma lower range and h having a 7 sigma upper range.
                                "--set-parameters", "cosmological_parameters", "omch2",     f"0.11336956221293837   {p['omegach2']}       0.12703091064106695",
                                "--set-parameters", "cosmological_parameters", "ombh2",     f"0.021615770756319073  {p['omegabh2']}     0.023103729722525095",
                                "--set-parameters", "cosmological_parameters", "h0",        f"0.6425290065540109    {p['h']}        0.7150188763839126",
                                "--set-parameters", "cosmological_parameters", "n_s",       f"0.9342806461291905    {p['ns']}       0.986698010466018",
                                "--set-parameters", "cosmological_parameters", "S_8_input", f"0.7229290272333152    {p['s8proxy']}     0.9134139168266714",
                                "--set-parameters", "cosmological_parameters", "tau",       f"0.015070795999054837  {p['tau']}     0.09381939280201967",
                                "--set-parameters", "planck",                  "a_planck",  f"0.9879083109867925    {p['calPlanck']}   1.0130810744845216",
                                "--set-priors", "planck", "a_planck", "gaussian 1.0 0.0025",]
                else:
                    cmd += ["--set-parameter", "cosmological_parameters", "omch2", f"nan {p['omegach2']} nan",
                            "--set-parameter", "cosmological_parameters", "ombh2", f"nan {p['omegabh2']} nan",
                            "--set-parameter", "cosmological_parameters", "h0", f"nan {p['h']} nan",]
                    if config_name not in ["fixAs", "fixAsns"]:
                        cmd += ["--set-parameter", "cosmological_parameters", "S_8_input", f"nan {p['s8proxy']} nan",]
                    if config_name not in ["fixns", "fixAsns"]:
                        cmd += ["--set-parameter", "cosmological_parameters", "n_s", f"nan {p['ns']} nan",]

                if "EE" in run_type:
                    cmd += ["--set-parameter", "intrinsic_alignment_parameters", "A", f"nan {p['a_ia']} nan",

                            "--set-parameter", "nofz_shifts", "p_1", f"nan {p['p_z1']} nan",
                            "--set-parameter", "nofz_shifts", "p_2", f"nan {p['p_z2']} nan",
                            "--set-parameter", "nofz_shifts", "p_3", f"nan {p['p_z3']} nan",
                            "--set-parameter", "nofz_shifts", "p_4", f"nan {p['p_z4']} nan",
                            "--set-parameter", "nofz_shifts", "p_5", f"nan {p['p_z5']} nan",]
                    if config_name in ["HMCode2020"]:
                        cmd += ["--set-parameters", "halo_model_parameters", "logT_AGN", f"7.3 {p['logt_agn']} 8.3",]
                    else:
                        cmd += ["--set-parameter", "halo_model_parameters", "A", f"nan {p['a_baryon']} nan",]

                if "nE" in run_type or "w" in run_type:
                    cmd += ["--set-parameter", "bias_parameters", "b1_bin_1", f"nan {p['b1l']} nan",
                            "--set-parameter", "bias_parameters", "b2_bin_1", f"nan {p['b2l']} nan",
                            "--set-parameter", "bias_parameters", "gamma3_bin_1", f"nan {p['gamma3l']} nan",
                            
    
                            "--set-parameter", "bias_parameters", "b1_bin_2", f"nan {p['b1h']} nan",
                            "--set-parameter", "bias_parameters", "b2_bin_2", f"nan {p['b2h']} nan",
                            "--set-parameter", "bias_parameters", "gamma3_bin_2", f"nan {p['gamma3h']} nan",]

                if "w" in run_type:
                    # avir is not used for GGL, so add it only if w is one of the probes.
                    cmd += ["--set-parameter", "bias_parameters", "a_vir_bin_1", f"nan {p['a_virl']} nan",
                            "--set-parameter", "bias_parameters", "a_vir_bin_2", f"nan {p['a_virh']} nan",]

                # dz prior means
                for i, m in enumerate(dx_mean):
                    cmd += ["--set-priors", "nofz_shifts", f"p_{i+1}", f"gaussian {m} 1.0"]

                # Config specific settings
                cmd += config
                if config_name == "wCDM":
                    cmd += ["--set-parameters", "cosmological_parameters", "w", f"-3.0 {p['w']} -0.33",]
                elif config_name == "oCDM":
                    cmd += ["--set-parameters", "cosmological_parameters", "omega_k", f"-0.4 {p['omegak']} 0.4",]
                elif config_name == "fRCDM":
                    cmd += ["--set-parameters", "cosmological_parameters", "log10_fR0", f"-8.0 {p['logfr0']} -2.0",]
                elif config_name == "nuCDM":
                    cmd += ["--set-parameters", "cosmological_parameters", "mnu", f"0.0 {p['mnu']} 3.0",]
                elif config_name == "Pantheon":
                    cmd += ["--set-parameters", "supernova_params",  "M",  f"-22.0  {p['sn_m']}   -17.0",]

                if "nE" in run_type:
                    cmd += nE_scale_cuts

                # Allow wedges to time out
                if "w" in run_type:
                    cmd += timeout_setttings

                # Hacky...
                if run_type == "EE_nE_w_Planck":
                    cmd += ["--run-type", "EE_nE_w"]
                elif run_type == "Planck":
                    cmd += ["--cut-modules", "wedges",
                            "--cut-modules", "BOSS_like",
                            "--fix-values", "bias_parameters"]
                    cmd += ["--run-type", "w"]
                else:
                    cmd += ["--run-type", run_type]


                cmd += ["--sampler-config", "max_iterations", "3000",
                        "--sampler-config", "maxlike_method", solver]

                if solver == "Powell":
                    if "Planck" in run_type:
                        cmd += ["--sampler-config", "step_size", "0.01",
                                "--sampler-config", "maxlike_tolerance", "0.00001"]
                    else:
                        cmd += ["--sampler-config", "step_size", "0.1",
                                "--sampler-config", "maxlike_tolerance", "0.0001"]
                elif solver == "L-BFGS-B":
                    MAP_covmat_file = os.path.join(root_dir, run_name, "chain", f"cov_estimate_{run_name}.txt")
                    cmd += ["--sampler-config", "grad_eps", "1e-3",
                            "--sampler-config", "gtol", "0.001",
                            "--sampler-config", "ftol", "1e-6",
                            "--sampler-config", "lbfgs_maxcor", "200",
                            "--sampler-config", "output_covmat", MAP_covmat_file]

                #cmd += ["--overwrite"]

                subprocess.run(["python", script] + cmd, check=True)
                    