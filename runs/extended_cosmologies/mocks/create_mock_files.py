import subprocess

if __name__ == "__main__":
    script = "utils/run_kcap.py"
    output_dir = "runs/extended_cosmologies/mocks/data/noisefree_fRCDM/"
    
    # Create f(R) mocks

    run_name_root = "fRCDM"
    run_type = "EE_fR"
    run_name = f"{run_name_root}_{run_type}"

    cmd = ["--create-mocks", "--noiseless-mocks",
            "--root-dir", output_dir,
            "--KiDS-data-file", "../K1000-data/Phase-1/twopoint/twoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad.fits",
            "--run-name", run_name,
            "--run-type", run_type,
            "--overwrite"]
    cmd += ["--set-parameters", "cosmological_parameters", "log10_fR0", "-5.0"]
    #cmd += ["--set-keys", "reaction", "massloop", "30"]

    subprocess.run(["python", script] + cmd, check=True)