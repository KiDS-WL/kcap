
import subprocess
import os

import numpy as np

import sys
sys.path.append("utils")
import process_chains


if __name__ == "__main__":
    script = "utils/run_kcap.py"

    multinest_chain_files = {"A" : {"EE_nE_w" : "runs/3x2pt/data_initial_cov_fast/cosmology/multinest_blindA_EE_nE_w/chain/samples_multinest_blindA_EE_nE_w.txt",
                                    "EE_nE"   : "runs/3x2pt/data_initial_cov_fast_w_timeout_bad_w_zs/cosmology/multinest_blindA_EE_nE/chain/samples_multinest_blindA_EE_nE.txt",},
                             "B" : {"EE_nE_w" : "runs/3x2pt/data_initial_cov_fast/cosmology/multinest_blindB_EE_nE_w/chain/samples_multinest_blindB_EE_nE_w.txt",
                                    "EE_nE"   : "runs/3x2pt/data_initial_cov_bad_w_zs/cosmology/multinest_blindB_EE_nE/chain/samples_multinest_blindB_EE_nE.txt",},
                             "C" : {"EE_nE_w" : "runs/3x2pt/data_initial_cov_fast/cosmology/multinest_blindC_EE_nE_w/chain/samples_multinest_blindC_EE_nE_w.txt",
                                    "EE_nE"   : "runs/3x2pt/data_initial_cov_bad_w_zs/cosmology/multinest_blindC_EE_nE/chain/samples_multinest_blindC_EE_nE.txt",},
                            }
    n_start_point = 12

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

    root_dir = "runs/3x2pt/data_initial_cov_MAP_L_BFGS_B_ftol6_eps3_multi/cosmology/"

    blinds = ["A", "B", "C"]
    run_types = ["EE_nE_w",]# "EE_nE", ]

    sampler = "maxlike"
    run_name_root = "MAP"

    solver = "L-BFGS-B"

    for blind in blinds:
        print(f"Blind {blind}")
        twopoint_file = twopoint_file_template.format(blind=blind)
        multinest_chain_file = multinest_chain_files[blind]

        for run_type in run_types:
            print(f"  Run type: {run_type}")

            # Get starting points
            chain = process_chains.load_chain(multinest_chain_file[run_type])

            param_names = [n.name for n in chain.getParamNames().names]

            d = chain.getParams()
            max_logpost_idx = (d.logpost).argmax()
            logpost_sort_idx = np.argsort(d.logpost)[::-1]
            
            start_points = []
            for idx in logpost_sort_idx[:n_start_point]:
                start_points.append({n : getattr(d, n)[idx] for n in param_names})

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
                print(f"Starting point {i} (logpost : {p['logpost']:.2f})")

                run_name = f"{run_name_root}_{i}_blind{blind}_{run_type}_{solver}"

                # np.savetxt(os.path.join(root_dir, run_name, "config", "powell_search_directions.txt"), search_directions)

                # Base setup
                cmd = [ "--root-dir", root_dir,
                        "--run-name", run_name,
                        "--run-type", run_type,
                        "--KiDS-data-file", twopoint_file,
                        "--dz-covariance-file", dz_cov_file,
                        "--BOSS-data-files", *boss_data_files,
                        "--BOSS-covariance-files", *boss_cov_files,
                        "--sampler", sampler,]
                        
                cmd += ["--set-parameter", "cosmological_parameters", "omch2", f"nan {p['omegach2']} nan",
                        "--set-parameter", "cosmological_parameters", "ombh2", f"nan {p['omegabh2']} nan",
                        "--set-parameter", "cosmological_parameters", "h0", f"nan {p['h']} nan",
                        "--set-parameter", "cosmological_parameters", "n_s", f"nan {p['ns']} nan",
                        "--set-parameter", "cosmological_parameters", "S_8_input", f"nan {p['s8']} nan",

                        "--set-parameter", "halo_model_parameters", "A", f"nan {p['a_baryon']} nan",
                        "--set-parameter", "intrinsic_alignment_parameters", "A", f"nan {p['a_ia']} nan",

                        "--set-parameter", "nofz_shifts", "p_1", f"nan {p['p_z1']} nan",
                        "--set-parameter", "nofz_shifts", "p_2", f"nan {p['p_z2']} nan",
                        "--set-parameter", "nofz_shifts", "p_3", f"nan {p['p_z3']} nan",
                        "--set-parameter", "nofz_shifts", "p_4", f"nan {p['p_z4']} nan",
                        "--set-parameter", "nofz_shifts", "p_5", f"nan {p['p_z5']} nan",

                        "--set-parameter", "bias_parameters", "b1_bin_1", f"nan {p['b1l']} nan",
                        "--set-parameter", "bias_parameters", "b2_bin_1", f"nan {p['b2l']} nan",
                        "--set-parameter", "bias_parameters", "gamma3_bin_1", f"nan {p['gamma3l']} nan",
                        

                        "--set-parameter", "bias_parameters", "b1_bin_2", f"nan {p['b1h']} nan",
                        "--set-parameter", "bias_parameters", "b2_bin_2", f"nan {p['b2h']} nan",
                        "--set-parameter", "bias_parameters", "gamma3_bin_2", f"nan {p['gamma3h']} nan",
                        ]

                if "w" in run_type:
                    # avir is not used for GGL, so add it only if w is one of the probes.
                    cmd += ["--set-parameter", "bias_parameters", "a_vir_bin_1", f"nan {p['a_virl']} nan",
                            "--set-parameter", "bias_parameters", "a_vir_bin_2", f"nan {p['a_virh']} nan",]

                # dz prior means
                for i, m in enumerate(dx_mean):
                    cmd += ["--set-priors", "nofz_shifts", f"p_{i+1}", f"gaussian {m} 1.0"]

                if "nE" in run_type:
                    cmd += nE_scale_cuts

                # Allow wedges to time out
                if "w" in run_type:
                    cmd += timeout_setttings

                cmd += ["--sampler-config", "max_iterations", "3000",
                        "--sampler-config", "maxlike_method", solver]

                if solver == "Powell":
                    cmd += ["--sampler-config", "step_size", "0.1",
                            "--sampler-config", "maxlike_tolerance", "0.0001"]
                elif solver == "L-BFGS-B":
                    MAP_covmat_file = os.path.join(root_dir, run_name, "chain", f"cov_estimate_{run_name}.txt")
                    cmd += ["--sampler-config", "grad_eps", "1e-3",
                            "--sampler-config", "gtol", "0.001",
                            "--sampler-config", "ftol", "1e-6",
                            "--sampler-config", "lbfgs_maxcor", "200",
                            "--sampler-config", "output_covmat", MAP_covmat_file]

                # cmd += ["--overwrite"]

                subprocess.run(["python", script] + cmd, check=True)
                    