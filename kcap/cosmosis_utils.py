import glob
import collections
import os
import re
import configparser
import io
import warnings
import copy
import shutil

import numpy as np

import cosmosis.runtime.config
import cosmosis.runtime.pipeline
import cosmosis.datablock

def dict_to_datablock(d={}):
    b = cosmosis.datablock.DataBlock()
    for section in d.keys():
        for name, value in d[section].items():
            b[section, name] = value

    return b

def create_pipeline(config, verbose=True):
    modules = []
    config_copy = copy.deepcopy(config)
    for module_name, c in config_copy.items():
        filename = c.pop("file")
        module = cosmosis.runtime.module.Module(module_name=module_name,
                                                file_path=filename)
        module.setup(dict_to_datablock({module_name : c}))
        modules.append(module)

    def pipeline(block):
        for m in modules:
            if verbose: print(f"Running {m.name}")
            status = m.execute(block)
            if status != 0:
                raise RuntimeError(f"Module {m.name} failed at execute.")
    
    return pipeline
    
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

def flatten_config(values, only_str_list_conversion=False):
    values = copy.deepcopy(values)
    for section, d in values.items():
        for key, value in d.items():
            if isinstance(value, (list, tuple)) and all(isinstance(l, str) for l in value):
                # Join list of strings
                d[key] = " ".join([v for v in value])
            if not only_str_list_conversion:
                if isinstance(value, (list, tuple)) and not all(isinstance(l, str) for l in value):
                    # Join other lists as well
                    d[key] = " ".join([str(v) for v in value])
                if value is True:
                    d[key] = "T"
                elif value is False:
                    d[key] = "F"
    return values


def emulate_configparser_interpolation(string, defaults, max_depth=100):
    """Emulate the configparser basic interpolation syntax.
    
    Interpolation happens by replacing occurrences of %(key)s in string by
    defaults[key].
    """
    def replacer(match):
        s = defaults[match.group("replace")]
        return s

    old_string = string
    for depth in range(max_depth):
        new_string = re.sub(pattern=r"%\((?P<replace>[a-z0-9_]+)\)s", 
                            repl=replacer,
                            string=old_string, flags=re.IGNORECASE)
        if new_string == old_string:
            break
        else:
            old_string = new_string

    return new_string


class CosmoSISPipelineFactory:
    def __init__(self, options, files=None):
        self._init_options = copy.deepcopy(options)
        
        self.base_config = self.create_base_config(**options)
        self.base_params = self.create_base_params()
        self.reset_config()
        self.reset_params()
        
        self.file_options_registry = self.base_config_data_files
            
    def reset_config(self):
        self.config = copy.deepcopy(self.base_config)
        
    def reset_params(self):
        self.params = copy.deepcopy(self.base_params)
        
    def update_config(self, config_update):
        for sec, sec_update in config_update.items():
            # Delete section if section update value is Nine
            if sec_update is None:
                del self.config[sec]
            else:
                for opt, val in copy.deepcopy(sec_update).items():
                    # Delete option if update value is None
                    if val is None and opt in self.config[sec]:
                        del self.config[sec][opt]
                        # Remove from the update dict
                        del sec_update[opt]
                # Else update the option with the new values
                self.config[sec].update(sec_update)
                
    def update_params(self, params_update):
        for k, u in params_update.items():
            self.params[k].update(u)
                
    def add_section(self, section, position):
        pass
    
    def run_pipeline(self, defaults=None):
        c = copy.deepcopy(self.config)
        p = copy.deepcopy(self.params)
        
        # Resolve %()s definitions with defaults dict
        for sec in c:
            for option, value in c[sec].items():
                if isinstance(value, str):
                    c[sec][option] = emulate_configparser_interpolation(value, defaults)

        # Discard range information if available
        for sec in p:
            for name, value in p[sec].items():
                if isinstance(value, list):
                    p[sec][name] = value[1]            

        pipeline = create_pipeline(c)
        block = dict_to_datablock(p)
        pipeline(block)

        return block

    
    def add_sampling_config(self, sampling_options):
        modules = sampling_options.pop("modules",
                                       self.config.keys())
        derived_parameters = sampling_options.pop("derived_parameters", 
                                                  self.derived_params())
        parameter_file = sampling_options.pop("parameter_file",
                                              "%(CONFIG_DIR)s/values.ini")
        prior_file = sampling_options.pop("prior_file", 
                                          "%(CONFIG_DIR)s/priors.ini")
                    
        self.sampling_config = self.create_base_sampling_config(
                                        modules=modules,
                                        derived_parameters=derived_parameters,
                                        parameter_file=parameter_file,
                                        prior_file=prior_file,
                                        **sampling_options)
        
        self.config = {**self.sampling_config, **self.config}
        
    def stage_files(self, root_dir, defaults, data_file_dirs=None, copy_data_files=False, create_multinest_dir=True):
        pipeline_config = copy.deepcopy(self.config)
        
        # Create directories
        config_dir = os.path.join(root_dir, "config")
        os.makedirs(config_dir, exist_ok=True)
        os.makedirs(os.path.join(root_dir, "output"), exist_ok=True)
        if create_multinest_dir:
            os.makedirs(os.path.join(root_dir, "output", "multinest"), exist_ok=True)
            
        if copy_data_files:
            if self.file_options_registry is None:
                raise ValueError("No files have been registered for staging.")
                
            data_file_dirs = data_file_dirs or {}

            data_dir = os.path.join(root_dir, "data")
            os.makedirs(data_dir, exist_ok=True)
            
            for sec, name in self.file_options_registry:
                filepaths = self.config[sec][name]
                # Resolve %()s in the paths
                filepaths = emulate_configparser_interpolation(filepaths, {**defaults, **data_file_dirs})
                # Allow for case of multiple files per option. E.g. n(z) files
                filepaths = filepaths.split(" ")
                new_filepaths = []
                
                staging_dir = os.path.join(data_dir, sec)
                os.makedirs(staging_dir, exist_ok=True)
                
                for filepath in filepaths:
                    filename = os.path.split(filepath)[1]
                    shutil.copy(filepath, os.path.join(staging_dir, filename))
                    new_filepaths.append("%(DATA_DIR)s/"+f"{sec}/{filename}")
                
                new_filepaths = " ".join(new_filepaths) if len(new_filepaths) > 1 else new_filepaths[0]
                
                pipeline_config[sec][name] = new_filepaths
                
        d = copy.deepcopy(defaults)
        if not any([name in d for name in ["root_dir", "ROOT_DIR"]]):
            d["ROOT_DIR"] = root_dir
        if not any([name in d for name in ["data_dir", "DATA_DIR"]]):
            d["DATA_DIR"] = "%(ROOT_DIR)s/data/"
        if not any([name in d for name in ["output_dir", "OUTPUT_DIR"]]):
            d["OUTPUT_DIR"] = "%(ROOT_DIR)s/output/"
        if not any([name in d for name in ["data_dir", "CONFIG_DIR"]]):
            d["CONFIG_DIR"] = "%(ROOT_DIR)s/config/"

        # Create ini files
        ini = configparser.ConfigParser()
        ini.read_dict(flatten_config({"DEFAULT" : d, **pipeline_config}))
        with open(os.path.join(config_dir, "pipeline.ini"), "w") as f:
            ini.write(f)

        ini = configparser.ConfigParser()
        ini.read_dict(flatten_config(self.params))
        with open(os.path.join(config_dir, "values.ini"), "w") as f:
            ini.write(f)

        ini = configparser.ConfigParser()
        ini.read_dict({})
        with open(os.path.join(config_dir, "priors.ini"), "w") as f:
            ini.write(f)
    
    def create_base_config(self, *args, **kwargs):
        raise NotImplementedError("create_base_config should be overwritten "
                                  "by a subclass.")
    
    @property
    def base_config_data_files(self):
        raise NotImplementedError("base_config_data_files should be overwritten "
                                  "by a subclass.")
        
    def create_base_params(self, *args, **kwargs):
        raise NotImplementedError("create_base_params should be overwritten "
                                  "by a subclass.")
        
    def create_base_sampling_config(self,
                                    modules,
                                    derived_parameters,
                                    parameter_file,
                                    prior_file,
                                    verbose,
                                    debug,
                                    sampler_name,
                                    run_name,
                                    max_iterations=10000,
                                    resume=True,
                                    live_points=250,
                                    nested_sampling_tolerance=0.1,
                                    multinest_efficiency=0.8,
                                    multinest_const_efficiency=False,
                                    emcee_walker=80,
                                    emcee_covariance_file="",
                                    maxlike_method="Nelder-Mead",
                                    maxlike_tolerance=1e-3,
                                    max_posterior=True,
                                    **extra_sampler_options,
                             ):
        config = {      "pipeline" :   {"modules"           : " ".join(modules),
                                        "values"            : parameter_file,
                                        "priors"            : prior_file,
                                        "likelihoods"       : "tsz_like",
                                        "extra_output"      : " ".join(derived_parameters),
                                        "quiet"             : "F" if verbose else "T",
                                        "timing"            : "T",
                                        "debug"             : "T" if debug else "F"},

                        "runtime"  :   {"sampler"           : sampler_name},

                        "output"   :   {"filename"          : "%(OUTPUT_DIR)s/samples_%(RUN_NAME)s.txt",
                                        "format"            : "text"},

                        "multinest" :  {"max_iterations"    : max_iterations,
                                        "multinest_outfile_root" : "%(OUTPUT_DIR)s/multinest/multinest_%(RUN_NAME)s_",
                                        "update_interval"   : 20,
                                        "resume"            : "T" if resume else "F",
                                        "live_points"       : live_points,
                                        "efficiency"        : multinest_efficiency,
                                        "tolerance"         : nested_sampling_tolerance,
                                        "constant_efficiency" : "T" if multinest_const_efficiency else "F"},

                        "emcee"      : {"walkers"           : emcee_walker,
                                        "samples"           : max_iterations,
                                        "covmat"            : emcee_covariance_file,
                                        "nsteps"            : 5},

                        "test" :       {"save_dir"          : "%(OUTPUT_DIR)s/data_block",
                                        "fatal_errors"      : "T",},

                        "maxlike" :    {"method"          : maxlike_method,
                                        "tolerance"       : maxlike_tolerance,
                                        "maxiter"         : max_iterations,
                                        "max_posterior"   : "T" if max_posterior else "F",
                                        "output_steps"    : "T",
                                        "flush_steps"     : 1,},
                        }
        return config
    
    def derived_params(self):
        derived_parameters=["cosmological_parameters/S_8",
                            "cosmological_parameters/sigma_8",
                            "cosmological_parameters/A_s",
                            "cosmological_parameters/omega_m",
                            "cosmological_parameters/omega_nu",
                            "cosmological_parameters/omega_lambda"]
        if "correlate_dz" in self.config:
            derived_parameters += self.config["correlate_dz"]["output_parameters"].split(" ")
            
        return derived_parameters


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
    def __init__(self, module_file="boltzmann/pycamb/camb_interface.py", base_path=r"%(CSL_PATH)s",
                       transfer_function=True, **kwargs):
        super().__init__(module_file, base_path)

        self.config.update({"mode"              : "transfer" if transfer_function else "background",
                            "lmax"              : "2500",
                            "feedback"          : "0",
                            "zmax"              : str(float(kwargs.pop("zmax", 6.0))),
                            "nz"                : str(int(kwargs.pop("nz", 100))),
                            "background_zmax"   : str(float(kwargs.pop("background_zmax", 6.0))),
                            "background_nz"     : str(int(kwargs.pop("background_nz", 600))),})
        self.config.update({k : str(v) for k,v in kwargs.items()})

        self.default_parameters = {}#"cosmological_parameters" : {"omega_m"  : "0.3",
        #                                                         "h0"       : "0.70",
        #                                                         "omega_b"  : "0.05",
        #                                                         "tau"      : "0.089",
        #                                                         "n_s"      : "0.96",
        #                                                         "omega_k"  : "0.0",
        #                                                         "w"       : "-1.0",
        #                                                         "wa"      : "0.0",
        #                                                         "A_s"      : "2.1e-9",}}

class HalofitModule(CosmoSISModule):
    module_name = "halofit"
    def __init__(self, module_file=r"boltzmann/halofit_takahashi/halofit_interface.so", base_path=r"%(CSL_PATH)s",
                       **kwargs):
        super().__init__(module_file, base_path)
        self.config.update({k : str(v) for k,v in kwargs.items()})

class HMCodeModule(CosmoSISModule):
    module_name = "HMCode"
    def __init__(self, module_file=r"structure/meadcb/mead_interface.so", base_path=r"%(CSL_PATH)s",
                       **kwargs):
        super().__init__(module_file, base_path)
        self.config.update({k : str(v) for k,v in kwargs.items()})     

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
        elif hm_mode.lower() == "hmx":
            self.config["ihm"] = "18"

        self.default_parameters = {"halo_model_parameters" : {"log10_Theat" : "7.8"}}
        

class LoadNofzModule(CosmoSISModule):
    module_name = "load_nofz"
    
    def __init__(self, nofz_file, output_section="shear",
                       module_file=r"number_density/load_nz/load_nz.py", base_path=r"%(CSL_PATH)s",):
        self.module_name = self.module_name + "_" + output_section
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

        translate = {"galaxy" : "position"}
        probes_translated = [(translate.get(p,p) for p in probe) for probe in probes]
        self.config.update({"-".join(probe_t) : "-".join(probe) for probe_t, probe in zip(probes_translated,probes)})

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
                       module_file=r"libcosebis_cl.so", base_path=r"%(COSEBIS_PATH)s"):
        super().__init__(module_file, base_path)

        probe = probes[0]
        if len(probes) > 1:
            warnings.warn(f"Only first probe ({probe}) will be used.")
        if not (probe[0].lower() == "shear" and probe[1].lower() == "shear"):
            raise ValueError("COSEBIS only support shear-shear.")

        self.config.update({"theta_min"            : str(float(theta_min)),
                            "theta_max"            : str(float(theta_max)),
                            "n_max"                : str(int(n_max)),
                            "input_section_name"   : self.input_section_name,
                            "output_section_name"  : self.output_section_name,
                            "is_it_bmodes"         : "1" if b_modes else "0",
                            "Wn_Output_FolderName" : r"%(COSEBIS_PATH)s/WnLog/",
                            "Roots_n_Norms_FolderName" : r"%(COSEBIS_PATH)s/TLogsRootsAndNorms",
                            "Tn_Output_FolderName"     : r"%(COSEBIS_PATH)s/TpnLog/"})

class BOSSModule(CosmoSISModule):
    module_name = "boss"
    output_section_name_wedges = "xi_wedges"
    
    def __init__(self, window_file="", bands_file="", z_eff=[0.61], verbose=False,
                       module_file=r"cosmosis_module.py", base_path=r"%(BOSS_PATH)s",
                       **kwargs):
        super().__init__(module_file, base_path)

        self.config.update({"window_file"          : window_file,
                            "bands_file"           : bands_file,
                            "z_eff"                : " ".join([str(z) for z in z_eff]),
                            "output_section_name_wedges"  : self.output_section_name_wedges,
                            "verbose"              : "T" if verbose else "F",
                            })
        self.config.update({k : str(v) for k,v in kwargs.items()})

class InterpolatePowerSpectrumModule(CosmoSISModule):
    module_name = "galaxy_bias"
    def __init__(self, module_file=r"utils/galaxy_bias.py", base_path=r"%(KCAP_PATH)s",
                       **kwargs):
        super().__init__(module_file, base_path)
        self.config.update({k : str(v) for k,v in kwargs.items()})

class CosmoSISPipeline:
    def __init__(self, paths={}, sampler="", parameter_file="", prior_file="", 
                       verbose=0, debug=False):
        self.paths = paths
        self.paths["KCAP_PATH"] = os.path.split(os.path.dirname(__file__))[0]
        if "CSL_PATH" not in paths:
            self.paths["CSL_PATH"] = os.path.join(self.paths["KCAP_PATH"], "cosmosis-standard-library")
        if "COSEBIS_PATH" not in paths:
            self.paths["COSEBIS_PATH"] = os.path.join(self.paths["KCAP_PATH"], "cosebis")
        if "BOSS_PATH" not in paths:
            self.paths["BOSS_PATH"] = os.path.join(self.paths["KCAP_PATH"], "boss")

        self.config_defaults = {key : path for key, path in self.paths.items()}

        self.cosmosis_config = self.setup_pipeline_config(sampler="",
                                   parameter_file=parameter_file, prior_file=prior_file, 
                                   verbose=verbose, debug=debug)

        self.modules = []

    def assemble_ini(self, parameters=None):
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
        
        return ini, values

    def create_pipeline(self, ini, values):
        cosmosis_ini = config_to_cosmosis_ini(ini)
        pipeline = cosmosis.runtime.pipeline.LikelihoodPipeline(cosmosis_ini, values=values)
        return pipeline

    def run(self, parameters={}):
        self.pipeline_ini, values = self.assemble_ini(parameters)
        self.pipeline = self.create_pipeline(self.pipeline_ini, values)
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