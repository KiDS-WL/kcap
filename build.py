import pkg_resources
import subprocess
import os
import sys
import warnings
import argparse

minimum_cosmosis_version = pkg_resources.parse_version("0.0.8")

minimum_cc_version = pkg_resources.parse_version("5.0.0")
minimum_cxx_version = pkg_resources.parse_version("5.0.0")

def check_cosmosis():
    try:
        import cosmosis
        try:
            cosmosis_version = cosmosis.__version__
        except AttributeError:
            raise ModuleNotFoundError(f"Installed CosmoSIS version does not meet requirements ({minimum_cosmosis_version}).")
        cosmosis_version = pkg_resources.parse_version(cosmosis_version)
        
        if minimum_cosmosis_version > cosmosis_version:
            raise ModuleNotFoundError(f"Installed CosmoSIS version ({cosmosis_version}) does not meet requirements ({minimum_cosmosis_version}).")
    except ModuleNotFoundError:
        env = check_compilers()
        install_cosmosis(env)

    import cosmosis
    cosmosis_src_dir = os.path.split(os.path.dirname(cosmosis.__file__))[0]
    return {"COSMOSIS_SRC_DIR" : cosmosis_src_dir}

def install_cosmosis(env):
    cosmosis_source = "git+https://bitbucket.org/tilmantroester/cosmosis.git@kcap#egg=cosmosis-standalone"
    subprocess.check_call([sys.executable, "-m", "pip", "install", cosmosis_source], env=env)

def check_compilers():
    default_cc    = "gcc"
    default_cxx   = "g++"
    default_fc    = "gfortran"
    default_mpifc = "" # Disable MPI be default, else set to "mpif90"

    env = {"PATH"  : os.environ["PATH"] if "PATH" in os.environ else "/usr/bin/",
           "CC"    : os.environ["CC"] if "CC" in os.environ else default_cc,
           "CXX"   : os.environ["CXX"] if "CXX" in os.environ else default_cxx,
           "FC"    : os.environ["FC"] if "FC" in os.environ else default_fc,
           "MPIFC" : os.environ["MPIFC"] if "MPIFC" in os.environ else default_mpifc,}

    cc_version = subprocess.check_output("${CC} -dumpversion", shell=True, env=env).decode("utf-8")
    cxx_version = subprocess.check_output("${CXX} -dumpversion", shell=True, env=env).decode("utf-8")
    fc_version = subprocess.check_output("${FC} -dumpversion", shell=True, env=env).decode("utf-8")

    cc_version = pkg_resources.parse_version(cc_version)
    cxx_version = pkg_resources.parse_version(cxx_version)

    if cc_version < minimum_cc_version:
        raise RuntimeError(f"GCC compiler version ({cc_version}) does not meet requirements ({minimum_cc_version}).")
    if cxx_version < minimum_cxx_version:
        raise RuntimeError(f"G++ compiler version ({cxx_version}) does not meet requirements ({minimum_cxx_version}).")
    
    if env["MPIFC"] == "":
        print("Compiling CosmoSIS without MPI support. If MPI support is required, set MPIFC.")
 
    return env

def build():
    # Check for sufficiently recent CosmoSIS-standalone installation
    cosmosis_env = check_cosmosis()

    os.makedirs("build", exist_ok=True)
    print("Running cmake.")
    subprocess.check_call(["cmake", ".."], cwd="build")
    print("Running make.")
    subprocess.check_call(["make"], cwd="build")
    print("Running make install.")
    subprocess.check_call(["make", "install"], cwd="build")

def clean():
    print("Running make csl-clean")
    subprocess.check_call(["make", "csl-clean"], cwd="build")
    print("Running make clean")
    subprocess.check_call(["make", "clean"], cwd="build")
    print("Removing CMakeCache.txt")
    os.remove("build/CMakeCache.txt")

def install_requirements(filename="requirements.txt"):
    print("Installing requirements.")
    with open(filename, "r") as f:
        lines = [l.rstrip() for l in f.readlines()]
        subprocess.check_call([sys.executable, "-m", "pip", "install"] + lines)

if __name__ == "__main__":
    # We're not using distutils/setuptools because we're not installing anything
    
    parser = argparse.ArgumentParser(description="Build the KiDS cosmology pipeline.", add_help=True)
    parser.add_argument("--clean", action='store_true', help="Clean up.")
    args = parser.parse_args(sys.argv[1:])

    if args.clean:
        clean()
    else:
        build()
        install_requirements()