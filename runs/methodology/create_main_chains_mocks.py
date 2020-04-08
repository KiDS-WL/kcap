import subprocess
import sys
import argparse
sys.path.append("modules/scale_cuts/")
from wrapper_twopoint import TwoPointWrapper

    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--noise-free', action='store_true', help='create noise-free mocks')
    parser.add_argument('--noise-range', nargs=2, default=[0, 0], type=int, metavar=('BEGIN', 'END'), help='create noisy mocks indexed by i with BEGIN <= i < END')
    parser.add_argument('--multi-start', nargs=2, default=[0, 0], type=int, metavar=('BEGIN', 'END'), help='create random staring point from multinest prior indexed by i with BEGIN <= i < END')
    args = parser.parse_args()

    script = "utils/run_kcap.py"
    
    data_name_root = "base"

    # Always create 3x2pt data vector for mocks
    # Configs know how to use only part of them properly
    run_type = "EE_nE_w"
    
    KiDS_twopoint_file = "../K1000-data/Phase-1/twopoint/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad.fits"
    
    # Always use identity matrix to generate mocks
    # Otherwise dz parameters will be transformed by the covariance
    dz_cov_file = "data/KV450/nofz/id_cov.asc"

    # Create noiseless mock data vector
    if args.noise_free:
        output_dir = "runs/methodology/data/noisefree_fiducial/"
        data_name = f"{data_name_root}_{run_type}"
        
        cmd = ["--create-mocks", "--noiseless-mocks",
                "--root-dir", output_dir,
                "--KiDS-data-file", KiDS_twopoint_file,
                "--dz-covariance-file", dz_cov_file,
                "--run-name", data_name,
                "--run-type", run_type]
        subprocess.run(["python", script] + cmd, check=True)

    # Create noisy mocks
    output_dir = "runs/methodology/data/noisy_fiducial/"

    for i in range(args.noise_range[0], args.noise_range[1]):
        data_name = f"{data_name_root}_{i}_{run_type}"

        cmd = ["--create-mocks",
                "--root-dir", output_dir,
                "--KiDS-data-file", KiDS_twopoint_file,
                "--dz-covariance-file", dz_cov_file,
                "--run-name", data_name,
                "--run-type", run_type]
        subprocess.run(["python", script] + cmd, check=True)


    # Multinest random
    multinest_prior_name = 'runs/methodology/data/multinest_prior/multinest_main_samples.npy'
    weight_name = 'runs/methodology/data/multinest_prior/multinest_main_weight.npy'
    data = np.load(multinest_prior_name)
    wgt  = np.load(weight_name)
    output_dir = f"runs/methodology/data/multinest_start/main_chain/"
    os.makedirs(f'{output_dir}', exist_ok=True)

    for i in range(args.multi_start[0], args.multi_start[1]):
        starting_point_settings = sample_random_from_multinest(data, wgt)
        np.save(f'{output_dir}start{i}.npy', starting_point_settings)

