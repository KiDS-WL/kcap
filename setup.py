import setuptools
import distutils.command.build

import pkg_resources
import subprocess
import os
import sys

minimum_cosmosis_version = pkg_resources.parse_version("0.0.8")

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


class CustomBuild(distutils.command.build.build):
    def run(self):
        # Check for sufficiently recent CosmoSIS-standalone installation
        cosmosis_env = check_cosmosis()

        super().run()

setuptools.setup(
    name="kcap",
    description="KiDS cosmology analysis pipeline",
    install_requires=["astropy"],
    cmdclass={"build" : CustomBuild},
)