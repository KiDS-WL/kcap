import subprocess

if __name__ == "__main__":

    # python utils/run_kcap.py --create-mocks --noiseless-mocks --root-dir data/KiDS1000/mocks/noisefree/ --KiDS-data-file runs/mocks/cosebis/COSEBIs_KiDS1000_omega_cam.fits --run-type cosebis --run-name cosebis_base --no-c-term --no-2d-c-term --use-nz-files --overwrite

    script = "utils/run_kcap.py"
    root_dir = "runs/mocks/cosebis/"
    
    base_twopoint_file = "data/KiDS1000/mocks/noisefree/cosebis_base/data/KiDS/COSEBIs_KiDS1000_omega_cam_mock_noiseless.fits"
    
    # Need to change
    base_dz_cov_file = "data/KV450/nofz/DIR_cov.asc"

    twopoint_file = base_twopoint_file
    dz_cov_file = base_dz_cov_file

    sampler = "multinest"
    
    run_type = "cosebis"
    run_name_root = "fast"
    run_name = f"{run_name_root}_{run_type}"
    cmd = ["--root-dir", root_dir,
           "--run-name", run_name,
           "--run-type", run_type,
           "--KiDS-data-file", twopoint_file,
           "--dz-covariance-file", dz_cov_file,
           "--sampler", sampler,
           "--sampler-config", "multinest_efficiency", "0.3",
           "--sampler-config", "nested_sampling_tolerance", "1.0e-2",
           "--cosebis-2d-c-term-file", "cosebis/example_files/cosebis/inputs/En_2D_cterm_KV450.ascii",
           "--cosebis-cos4phi-file", "cosebis/example_files/cosebis/inputs/En_cos4phi_KV450.ascii",
           "--cosebis-sin4phi-file", "cosebis/example_files/cosebis/inputs/En_sin4phi_KV450.ascii",
           "--use-nz-files",
           "--overwrite"]
    subprocess.run(["python", script] + cmd, check=True)