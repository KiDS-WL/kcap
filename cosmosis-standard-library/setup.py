import setuptools
import pkg_resources
import distutils.command.build
import distutils.command.install
import distutils.command.clean

import ctypes.util

import os
import subprocess
import sys
import warnings


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
        install_cosmosis()

    import cosmosis
    cosmosis_src_dir = os.path.split(os.path.dirname(cosmosis.__file__))[0]
    return {"COSMOSIS_SRC_DIR" : cosmosis_src_dir}

def install_cosmosis():
    cosmosis_source = "git+https://bitbucket.org/tilmantroester/cosmosis.git@standalone_setup_fixes#egg=cosmosis-standalone"
    subprocess.check_call([sys.executable, "-m", "pip", "install", cosmosis_source])

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

def find_library(lib, relative_include_path="include"):
    lib_path = ctypes.util.find_library(lib)
    inc_path = None
    if lib_path is not None:
        lib_path = os.path.dirname(lib_path)
        if os.path.split(lib_path)[1] == "lib":
            inc_path = os.path.join(os.path.split(lib_path)[0], relative_include_path)
        return lib_path, inc_path
    else:
        warnings.warn("Could not find {}.".format(lib))
        return "", ""
    
    
    

def check_libraries():
    env = {}
    # GSL
    if "GSL_LIB" not in os.environ and "GSL_INC" not in os.environ:
        gsl_lib_path, gsl_inc_path = find_library("gsl", "include/gsl")
        if gsl_lib_path is None or gsl_inc_path is None:
            print("GSL_LIB and GSL_INC not set. The GSL library and header directories can't be found on the system either.")
            print("Trying to install with conda.")
            subprocess.check_call(["conda", "install", "gsl"])
            gsl_lib_path, gsl_inc_path = find_library("gsl", "include/gsl")
            if gsl_lib_path is None or gsl_inc_path is None:
                raise RuntimeError("Still can't find GSL. Aborting")

        env["GSL_LIB"] = gsl_lib_path
        env["GSL_INC"] = gsl_inc_path
    else:
        env["GSL_LIB"] = os.environ["GSL_LIB"]
        env["GSL_INC"] = os.environ["GSL_INC"]

    # CFITSIO
    if "CFITSIO_LIB" not in os.environ and "CFITSIO_INC" not in os.environ:
        cfitsio_lib_path, cfitsio_inc_path = find_library("cfitsio", "include")
        env["CFITSIO_LIB"] = cfitsio_lib_path
        env["CFITSIO_INC"] = cfitsio_inc_path
    else:
        env["CFITSIO_LIB"] = os.environ["CFITSIO_LIB"]
        env["CFITSIO_INC"] = os.environ["CFITSIO_INC"]

    # FFTW
    if "FFTW_LIBRARY" not in os.environ and "FFTW_INCLUDE_DIR" not in os.environ:
        fftw_lib_path, fftw_inc_path = find_library("fftw3", "include/fftw3")
        env["FFTW_LIBRARY"] = fftw_lib_path
        env["FFTW_INCLUDE_DIR"] = fftw_inc_path
    else:
        env["FFTW_LIBRARY"] = os.environ["FFTW_LIBRARY"]
        env["FFTW_INCLUDE_DIR"] = os.environ["FFTW_INCLUDE_DIR"]

    # LAPACK
    if "LAPACK_LINK" not in os.environ:
        if sys.platform == "darwin":
            env["LAPACK_LINK"] = "-framework Accelerate"
        else:
            raise NotImplementedError(f"Platform {sys.platform} not supported yet.")
    else:
        env["LAPACK_LINK"] = os.environ["LAPACK_LINK"]

    # GFortran runtime library
    if "GFORTRAN_LIBPATH" not in os.environ:
        gfortran_lib_path, _ = find_library("gfortran")
        env["GFORTRAN_LIBPATH"] = gfortran_lib_path
    else:
        env["GFORTRAN_LIBPATH"] = os.environ["GFORTRAN_LIBPATH"]
    
    return env
            




class CustomBuild(distutils.command.build.build):
    def run(self):
        # Check for sufficiently recent CosmoSIS-standalone installation
        cosmosis_env = check_cosmosis()
        compiler_env = check_compilers()
        lib_env = check_libraries()

        env = {**cosmosis_env, **compiler_env, **lib_env}

        subprocess.check_call(["make"], env=env, cwd=".")

        super().run()

class CustomClean(distutils.command.clean.clean):
    def run(self):
        cosmosis_env = check_cosmosis()
        lib_env = check_libraries()

        env = {**cosmosis_env, **lib_env}
        subprocess.check_call(["make", "clean"], env=env, cwd=".")
        super().run()

class CustomInstall(distutils.command.install.install):
    def run(self):
        print("You probably don't want to install this...")


setuptools.setup(
    name="cosmosis-standalone-csl",
    description="CosmoSIS standalone standard library",
    install_requires=["cosmosis-standalone>=0.0.8"],
    cmdclass={"build" : CustomBuild,
              "build_ext" : CustomBuild,
              "clean" : CustomClean,
              "install" : CustomInstall,},
)