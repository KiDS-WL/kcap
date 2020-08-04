
import subprocess
import os

import numpy as np

import sys
sys.path.append("utils")
import process_chains


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

    root_dir = "runs/3x2pt/data_iterated_cov_MAP/test_sampler_MAP_run/"

    blinds = ["C"]
    run_types = ["EE_nE_w",]# "EE_nE", ]

    sampler = "test"
    run_name_root = "test"

    for blind in blinds:
        print(f"Blind {blind}")
        twopoint_file = twopoint_file_template.format(blind=blind)

        for run_type in run_types:
            print(f"  Run type: {run_type}")
            run_name = f"{run_name_root}_blind{blind}_{run_type}"

            # Base setup
            cmd = [ "--root-dir", root_dir,
                    "--run-name", run_name,
                    "--run-type", run_type,
                    "--KiDS-data-file", twopoint_file,
                    "--dz-covariance-file", dz_cov_file,
                    "--BOSS-data-files", *boss_data_files,
                    "--BOSS-covariance-files", *boss_cov_files,
                    "--sampler", sampler,]

            # dz prior means
            for i, m in enumerate(dx_mean):
                cmd += ["--set-priors", "nofz_shifts", f"p_{i+1}", f"gaussian {m} 1.0"]

            if "nE" in run_type:
                cmd += nE_scale_cuts

            # Allow wedges to time out
            if "w" in run_type:
                cmd += timeout_setttings

            # cmd += ["--overwrite"]

            subprocess.run(["python", script] + cmd, check=True)
                