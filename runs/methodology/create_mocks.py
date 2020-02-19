import sys
import subprocess

if __name__ == "__main__":

    if len(sys.argv) == 1:
        raise AssertionError('Usage:\n  $ python create_mocks.py main_chains\n  $ python create_mocks.py test_covariance')
    
    if len(sys.argv) >= 4:
        noise_begin = int(sys.argv[2])
        noise_end   = int(sys.argv[3])
    else:
        noise_begin = 0
        noise_end   = 1

    script = "utils/run_kcap.py"
    run_type = "EE_nE_w" # Create 3x2pt data vector
    
    if sys.argv[1] == 'main_chains':

        KiDS_twopoint_file =  "../K1000-data/Phase-1/twopoint/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad.fits"
        run_name_root = "base"

        # Create noiseless mock data vector
        if noise_begin == 0:
            output_dir = "runs/methodology/data/noisefree_fiducial/"
            run_name = f"{run_name_root}_{run_type}"
            cmd = ["--create-mocks", "--noiseless-mocks",
                    "--root-dir", output_dir,
                    "--KiDS-data-file", KiDS_twopoint_file,
                    "--run-name", run_name,
                    "--run-type", run_type]
            subprocess.run(["python", script] + cmd, check=True)

        # Create noisy mocks
        output_dir = "runs/methodology/data/noisy_fiducial/"

        for i in range(noise_begin, noise_end):
            run_name = f"{run_name_root}_{i}_{run_type}"

            cmd = ["--create-mocks",
                    "--root-dir", output_dir,
                    "--KiDS-data-file", KiDS_twopoint_file,
                    "--run-name", run_name,
                    "--run-type", run_type]
            subprocess.run(["python", script] + cmd, check=True)

    elif sys.argv[1] == 'test_covariance':

        print('WARNING: you might need to have run main_chains to avoid errors')
        
        KiDS_twopoint_tag_list = ['theoryEgretta']
        run_name_root_list     = ['theory_complex']

        # Create noiseless mock data vector
        for KiDS_twopoint_tag, run_name_root in zip(KiDS_twopoint_tag_list, run_name_root_list):
            KiDS_twopoint_file =  f"../K1000-data/Phase-1/twopoint/twoPoint_PneE+PeeE_mean_None_cov_{KiDS_twopoint_tag}_nOfZ_bucerosBroad.fits"

            # Create noiseless mock data vector
            if noise_begin == 0:
                output_dir = "runs/methodology/data/noisefree_fiducial/"
                run_name = f"{run_name_root}_{run_type}"

                cmd = ["--create-mocks", "--noiseless-mocks",
                        "--root-dir", output_dir,
                        "--KiDS-data-file", KiDS_twopoint_file,
                        "--run-name", run_name,
                        "--run-type", run_type]
                subprocess.run(["python", script] + cmd, check=True)

            # Create noisy mocks
            # We want to use the same noisy mocks as `main_chains`
            # So we run `wrapper_twopoint.py` instead of `run_kcap.py` here
            script2 = "modules/scale_cuts/wrapper_twopoint.py"
            output_dir = "runs/methodology/data/noisy_fiducial/"
            KiDS_twopoint_file_cov = KiDS_twopoint_file

            for i in range(noise_begin, noise_end):
                run_name = f"{run_name_root}_{i}_{run_type}"

                KiDS_twopoint_file_mean = f"{output_dir}base_{i}_EE_nE_w/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits"
                KiDS_twopoint_file_save = f"{output_dir}{run_name}/data/KiDS/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad_mock_noisy.fits"
                
                subprocess.run(['mkdir', '-p', f'{output_dir}{run_name}/data/KiDS/'])
                cmd = ['new_from_mean_and_cov', KiDS_twopoint_file_mean, KiDS_twopoint_file_cov, KiDS_twopoint_file_save]
                subprocess.run(["python", script2] + cmd, check=True)

