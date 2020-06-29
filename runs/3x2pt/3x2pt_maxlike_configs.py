
import subprocess
import os

import numpy as np

import sys
sys.path.append("utils")
import process_chains


if __name__ == "__main__":
    script = "utils/run_kcap.py"

    multinest_chain_files = {"A" : {"EE_nE_w" : "runs/3x2pt/data_initial_cov_fast_w_timeout/cosmology/multinest_blindA_EE_nE_w/chain/samples_multinest_blindA_EE_nE_w.txt"},
                             "B" : {"EE_nE_w" : "runs/3x2pt/data_initial_cov/cosmology/multinest_blindB_EE_nE_w/chain/samples_multinest_blindB_EE_nE_w.txt"},
                             "C" : {"EE_nE_w" : "runs/3x2pt/data_initial_cov/cosmology/multinest_blindC_EE_nE_w/chain/samples_multinest_blindC_EE_nE_w.txt"},
                            }
    n_start_point = 24

    # Cosmic shear + GGL twopoint file
    twopoint_file_template = "../Cat_to_Obs_K1000_P1/data/kids/fits/bp_KIDS1000_Blind{blind}_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits"

    # Covariance of the n(z)
    dz_cov_file = "../Cat_to_Obs_K1000_P1/data/kids/nofz/SOM_cov_multiplied.asc"

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

    root_dir = "runs/3x2pt/data_initial_cov_MAP/cosmology/"

    blinds = ["A", "B", "C"]
    run_types = ["EE_nE_w",]

    sampler = "maxlike"
    run_name_root = "MAP"

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

            for i, p in enumerate(start_points):
                print(f"Starting point {i} (logpost : {p['logpost']:.2f})")

                run_name = f"{run_name_root}_{i}_blind{blind}_{run_type}"

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
                        "--set-parameter", "bias_parameters", "a_vir_bin_1", f"nan {p['a_virl']} nan",

                        "--set-parameter", "bias_parameters", "b1_bin_2", f"nan {p['b1h']} nan",
                        "--set-parameter", "bias_parameters", "b2_bin_2", f"nan {p['b2h']} nan",
                        "--set-parameter", "bias_parameters", "gamma3_bin_2", f"nan {p['gamma3h']} nan",
                        "--set-parameter", "bias_parameters", "a_vir_bin_2", f"nan {p['a_virh']} nan",]

                if "nE" in run_type:
                    cmd += nE_scale_cuts

                # Allow wedges to time out
                if "w" in run_type:
                    cmd += timeout_setttings

                cmd += ["--sampler-config", "max_iterations", "3000",]

                # cmd += ["--overwrite"]

                subprocess.run(["python", script] + cmd, check=True)
                    