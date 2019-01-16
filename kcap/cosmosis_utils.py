import glob
import collections
import os
import configparser
import io
import warnings

import numpy as np

import cosmosis.runtime.config
import cosmosis.runtime.pipeline

def config_to_cosmosis_ini(config):
    ini = cosmosis.runtime.config.Inifile(None, )
    with io.StringIO() as s:
        config.write(s)
        s.seek(0)
        ini.read_file(s)
    return ini

def config_to_string(config):
    with io.StringIO() as s:
        config.write(s)
        s.seek(0)
        return s.read()

class CosmoSISModule:
    module_name = None
    def __init__(self, module_file, base_path=None):
        self.module_file = module_file
        if base_path:
            self.module_file =  os.path.join(base_path, module_file)

        self.config = collections.OrderedDict()
        self.config["file"] = self.module_file

        self.default_parameters = {}

    def get_section(self):
        return {self.module_name : self.config}

class ConsistencyModule(CosmoSISModule):
    module_name = "consistency"
    def __init__(self, module_file="utility/consistency/consistency_interface.py", base_path=r"%(CSL_PATH)s"):
        super().__init__(module_file, base_path)

class CAMBModule(CosmoSISModule):
    module_name = "camb"
    def __init__(self, module_file="boltzmann/camb/camb.so", base_path=r"%(CSL_PATH)s",
                       z_min=0.0, z_max=6.0, n_z=100, transfer_function=True):
        super().__init__(module_file, base_path)

        self.config.update({"mode"              : "all" if transfer_function else "background",
                            "lmax"              : "2500",
                            "feedback"          : "0",
                            "zmax"              : str(float(z_max)),
                            "nz"                : str(n_z),
                            "background_zmax"   : str(float(z_max)),
                            "background_nz"     : str(int(z_max*100)),})

        self.default_parameters = {"cosmological_parameters" : {"omega_m"  : "0.3",
                                                                "h0"       : "0.70",
                                                                "omega_b"  : "0.05",
                                                                "tau"      : "0.089",
                                                                "n_s"      : "0.96",
                                                                "omega_k"  : "0.0",
                                                                "w"       : "-1.0",
                                                                "wa"      : "0.0",
                                                                "A_s"      : "2.1e-9",}}

class HMxModule(CosmoSISModule):
    module_name = "HMx"
    def __init__(self, module_file=r"lib/cosmosis_interface.so", base_path=r"%(HMX_PATH)s",
                       transfer_function="camb",
                       k_min=1e-3, k_max=1e1, n_k=128, 
                       z_min=0.0, z_max=5.999, n_z=32,
                       fields=["dmonly"],
                       matter_matter_section_name="matter_power_nl",
                       verbose=0, hm_mode="hmcode", one_parameter_hmcode=False,
                       ihm=None):
        super().__init__(module_file, base_path)

        self.config.update({"fields"            : " ".join(fields),
                            "nk"                : str(n_k),
                            "kmin"              : str(float(k_min)),
                            "kmax"              : str(float(k_max)),
                            "nz"                : str(n_z),
                            "zmin"              : str(float(z_min)),
                            "zmax"              : str(float(z_max)),
                            "verbose"           : str(verbose),
                            "dimensionless_power_spectrum" : 0,
                            "p_lin_source"      : "external" if transfer_function.lower() == "camb" else "eh",
                            "matter_matter_section_name" : matter_matter_section_name,
                            "hm_mode"           : hm_mode,
                            "one_parameter_hmcode" : "T" if one_parameter_hmcode else "F"})
        if ihm is not None:
            self.config["ihm"] = str(ihm)

        self.default_parameters = {"halo_model_parameters" : {"log10_Theat" : "7.8"}}
        
class LoadNofzModule(CosmoSISModule):
    module_name = "load_nofz"
    
    def __init__(self, nofz_file, output_section="shear",
                       module_file=r"number_density/load_nz/load_nz.py", base_path=r"%(CSL_PATH)s",):
        super().__init__(module_file, base_path)
        
        if isinstance(nofz_file, str):
            filepath = nofz_file
        elif isinstance(nofz_file, (list, tuple)):
            filepath = " ".join(nofz_file)
        else:
            raise ValueError("nofz_file needs to be str or list of str.")
        self.config.update({"filepath"          : filepath,
                            "histogram"         : "True",
                            "output_section"    : f"nz_{output_section}",})

class ProjectionModule(CosmoSISModule):
    module_name = "projection"
    
    def __init__(self, probes=[], n_ell=400, ell_min=0.1, ell_max=5e4, verbose=True,
                       module_file=r"structure/projection/project_2d.py", base_path=r"%(CSL_PATH)s"):
        super().__init__(module_file, base_path)
     
        self.config.update({"ell_min"           : str(float(ell_min)),
                            "ell_max"           : str(float(ell_max)),
                            "n_ell"             : str(n_ell),
                            "verbose"           : "T" if verbose else "F",
                            "get_kernel_peaks"  : "F"})
        self.config.update({"-".join(probe) : "-".join(probe) for probe in probes})

class Cl2xiModule(CosmoSISModule):
    module_name = "cl2xi"
    
    def __init__(self, probes=[], 
                       module_file=r"shear/cl_to_xi_nicaea/nicaea_interface.so", base_path=r"%(CSL_PATH)s"):
        super().__init__(module_file, base_path)

        probe = probes[0]
        if len(probes) > 1:
            warnings.warn(f"Only first probe ({probe}) will be used.")
        if probe[0].lower() == "shear" and probe[1].lower() == "shear":
            corr_type = 0
        elif probe[0].lower() == "shear" or probe[1].lower() == "shear":
            corr_type = 1
        else:
            corr_type = 2
        self.config.update({"corr_type"         : str(corr_type)})

class COSEBISModule(CosmoSISModule):
    module_name = "cosebis"
    input_section_name = "shear_cl"
    output_section_name = "cosebis"
    
    def __init__(self, probes=[], theta_min=0.5, theta_max=300.0, n_max=5, b_modes=False,
                       module_file=r"cosebis/libcosebis_cl.so", base_path=r"%(KCAP_PATH)s"):
        super().__init__(module_file, base_path)

        probe = probes[0]
        if len(probes) > 1:
            warnings.warn(f"Only first probe ({probe}) will be used.")
        if not (probe[0].lower() == "shear" and probe[1].lower() == "shear"):
            raise ValueError("COSEBIS only support shear-shear.")

        self.config.update({"theta_min"           : str(float(theta_min)),
                            "theta_max"           : str(float(theta_max)),
                            "n_max"               : str(int(n_max)),
                            "input_section_name" : self.input_section_name,
                            "output_section_name" : self.output_section_name,
                            "is_it_bmodes"        : "1" if b_modes else "0",})

class CosmoSISPipeline:
    def __init__(self, paths={}, sampler="", parameter_file="", prior_file="", 
                       verbose=0, debug=False):
        self.paths = paths
        self.paths["KCAP_PATH"] = os.path.split(os.path.dirname(__file__))[0]
        if "CSL_PATH" not in paths:
            self.paths["CSL_PATH"] = os.path.join(self.paths["KCAP_PATH"], "cosmosis-standard-library")

        self.config_defaults = {key : path for key, path in self.paths.items()}

        self.cosmosis_config = self.setup_pipeline_config(sampler="",
                                   parameter_file=parameter_file, prior_file=prior_file, 
                                   verbose=verbose, debug=debug)

        self.modules = []

    def create_pipeline(self, parameters=None):
        ini = configparser.ConfigParser(defaults=self.config_defaults)
        ini.read_dict(self.cosmosis_config)

        ini["pipeline"]["modules"] = " ".join(m.module_name for m in self.modules)
        for m in self.modules:
            ini[m.module_name] = m.config

        if ini["pipeline"]["values"] == "":
            values = configparser.ConfigParser()
            for m in self.modules:
                values.read_dict(m.default_parameters)
            if parameters is not None:
                values.read_dict(parameters)
            values = config_to_cosmosis_ini(values)
        else:
            values = None
        
        # print("INI", config_to_string(ini))
        cosmosis_ini = config_to_cosmosis_ini(ini)
        pipeline = cosmosis.runtime.pipeline.LikelihoodPipeline(cosmosis_ini, values=values)
        return pipeline, ini

    def run(self, parameters={}):
        self.pipeline, self.pipeline_ini = self.create_pipeline()
        p = [p.start for p in self.pipeline.parameters]
        for section in parameters.keys():
            for name, value in parameters[section].items():
                idx = self.pipeline.parameter_index(section.lower(), name.lower())
                p[idx] = value
        data = self.pipeline.run_parameters(p, all_params=True)
        return data

    # def setup_base_config(self, output_path=None,
    #                       sampler="", sampler_output_file=None, 
    #                       project_triad_path=None, **sampler_kwargs):
    #     cosmosis_config = configparser.ConfigParser()  

    #     cosmosis_config["runtime"] =           {"sampler"               : sampler}

    #     # if sampler_output_file is None:
    #     #     sampler_output_file = os.path.join(output_path, "sampler_output.txt")

    #     # cosmosis_config["output"] =            {"format"                : "text",
    #     #                                         "filename"              : sampler_output_file}
    #     # if sampler.lower() == "test":
    #     #     cosmosis_config["test"] =          {"save_dir"               : r"%(OUTPUT_PATH)s",
    #     #                                         "fatal_errors"           : "T"}
    #     # if sampler.lower() == "fisher":
    #     #     cosmosis_config["fisher"] = sampler_kwargs

    #     return cosmosis_config

    def setup_pipeline_config(self, sampler="",
                                    parameter_file="", prior_file="", 
                                    verbose=0, debug=False):
        cosmosis_config = configparser.ConfigParser()  

        cosmosis_config["runtime"] =     {"sampler"           : sampler}

        cosmosis_config["pipeline"] =    {"values"            : parameter_file,
                                          "priors"            : prior_file,
                                          "likelihoods"       : "",
                                          "extra_output"      : "",
                                          "quiet"             : "T" if verbose == 0 else "F",
                                          "timing"            : "T" if verbose > 1 else "F",
                                          "debug"             : "T" if debug else "F"}
        return cosmosis_config


def load_cosmosis_params(output_path, section="cosmological_parameters"):
    with open(os.path.join(output_path, section, "values.txt"), "r") as f:
        params = {line.split("=")[0].strip() : float(line.split("=")[1].strip()) for line in f.readlines()}
    return params

def load_cosmosis_2pt(output_path, x_name="ell", y_name="shear", suffix="cl"):
    cosmosis_2pt = {}

    name = y_name + "_" + suffix
    for path in glob.glob(os.path.join(output_path, f"{name}*")):
        section = os.path.split(path)[1]
        cosmosis_2pt[section] = {}

        values = load_cosmosis_params(os.path.split(path)[0], section=section)
        cosmosis_2pt[section] = {**values}
        
        if os.path.isfile(os.path.join(path, x_name+".txt")):
            x_file = True
            x = np.loadtxt(os.path.join(path, x_name+".txt"))
        else:
            x_file = False
            
        for bin_path in glob.glob(os.path.join(path, "*_*_*")):
            key = os.path.split(bin_path)[1]
            key = key[:key.rfind(".txt")]
            data = np.loadtxt(bin_path).T
            if x_file:
                y = data
            else:
                x = data[0]
                y = data[1]
            
            cosmosis_2pt[section][key] = y
            
        cosmosis_2pt[section][x_name] = x
        
    return cosmosis_2pt