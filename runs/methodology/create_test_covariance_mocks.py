import subprocess
import os
import argparse
import sys
sys.path.append("modules/scale_cuts/")
from wrapper_twopoint import TwoPointWrapper

    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--noise-free', action='store_true', help='create noise-free mocks')
    parser.add_argument('--noise-range', nargs=2, default=[0, 1], type=int, metavar=('BEGIN', 'END'), help='create noisy mocks indexed by i with BEGIN <= i < END')
    args = parser.parse_args()

    if args.noise_range[1] > args.noise_range[0]:
        print('You want to create noisy mocks.')
        print('Note that noisy mocks for test_covariance require those for main_chains.')
        print('Make sure you have already created them.')

    script = "utils/run_kcap.py"

    # Always create 3x2pt data vector for mocks
    # Configs know how to use only part of them properly
    run_type = "EE_nE_w"

    KiDS_twopoint_tag_list = ['theoryBuceros', 'theoryEgretta', 'simBuceros', 'simEgretta']
    data_name_root_list    = ['theory_simple', 'theory_complex', 'mock_simple', 'mock_complex']

    # Create noiseless mock data vector
    for KiDS_twopoint_tag, data_name_root in zip(KiDS_twopoint_tag_list, data_name_root_list):
        KiDS_twopoint_file =  f"../K1000-data/Phase-1/twopoint/twoPoint_PneE+PeeE_mean_None_cov_{KiDS_twopoint_tag}_nOfZ_bucerosBroad.fits"

        # Create noiseless mock data vector
        if args.noise_free:
            output_dir = "runs/methodology/data/noisefree_fiducial/"
            data_name = f"{data_name_root}_{run_type}"

            cmd = ["--create-mocks", "--noiseless-mocks",
                    "--root-dir", output_dir,
                    "--KiDS-data-file", KiDS_twopoint_file,
                    "--run-name", data_name,
                    "--run-type", run_type]
            subprocess.run(["python", script] + cmd, check=True)

        # Create noisy mocks
        # We want to use the same noisy mocks as `main_chains`
        # So we create new files by combining old mocks with new covariance
        output_dir = "runs/methodology/data/noisy_fiducial/"
        TP_cov = TwoPointWrapper.from_fits(KiDS_twopoint_file, covmat_name='COVMAT')

        for i in range(args.noise_range[0], args.noise_range[1]):
            data_name = f"{data_name_root}_{i}_{run_type}"

            KiDS_twopoint_file_mean = f"{output_dir}base_{i}_EE_nE_w/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits"
            KiDS_twopoint_file_save = f"{output_dir}{data_name}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_{KiDS_twopoint_tag}_nOfZ_bucerosBroad_mock_noisy.fits"

            os.makedirs(f'{output_dir}{data_name}/data/KiDS/', exists_ok=True)
            TP_mean = TwoPointWrapper.from_fits(KiDS_twopoint_file_mean)
            mean    = TP_mean.makeMeanVector()
            TP_cov.replaceMeanVector(mean)
            TP_cov.to_fits(KiDS_twopoint_file_save, overwrite=True)

