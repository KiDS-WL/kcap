import setuptools
import distutils.command.build
import distutils.command.clean
import distutils.command.install

import pkg_resources
import subprocess
import os
import sys

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
    if "MPIFC" in env:
        print("Compiling with MPI support.")
    else:
        print("Compiling without MPI support.")
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
    
    return env

class CustomBuild(distutils.command.build.build):
    def run(self):
        # Check for sufficiently recent CosmoSIS-standalone installation
        cosmosis_env = check_cosmosis()

        os.makedirs("build", exist_ok=True)
        print("Running cmake.")
        subprocess.check_call(["cmake", ".."], cwd="build")
        print("Running make.")
        subprocess.check_call(["make"], cwd="build")
        print("Running make install.")
        subprocess.check_call(["make", "install"], cwd="build")

        super().run()

class CustomClean(distutils.command.clean.clean):
    def run(self):
        print("Running make csl-clean")
        subprocess.check_call(["make", "csl-clean"], cwd="build")
        print("Running make clean")
        subprocess.check_call(["make", "clean"], cwd="build")
        print("Removing CMakeCache.txt")
        os.remove("build/CMakeCache.txt")

        super().run()

class CustomInstall(distutils.command.install.install):
    def run(self):
        print("You probably don't want to install this...")

setuptools.setup(
    name="kcap",
    description="KiDS cosmology analysis pipeline",
    install_requires=["astropy"],
    cmdclass={"build"   : CustomBuild,
              "install" : CustomInstall,
              "clean"   : CustomClean,
              "develop" : CustomInstall},
)