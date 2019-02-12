import collections
import os
import configparser
import io

import numpy as np
from .cosmosis_utils import CosmoSISPipeline, ConsistencyModule, CAMBModule, \
                            HMxModule, ProjectionModule, LoadNofzModule, \
                            Cl2xiModule, COSEBISModule, \
                            config_to_string

def dict_insert(d, idx, obj, after=False):
    if after: idx += 1
    d_list = list(d.items())
    return collections.OrderedDict(d_list[:idx] + [obj] + d_list[idx:])

class TwoPoint:
    supported_fields_2D = ("shear", "y")
    supported_fields_3D = ("matter", "pressure")

    def __init__(self, probes={}, parameters={},
                 transfer_function="camb", nonlinear_Pk="HMx", 
                 tomographic_bin_names={}, nofz_files={},
                 module_configs={},
                 module_paths={}):

        self.probes = probes

        statistics = []
        fields_3D = []
        fields_2D = []
        for statistic, fields in probes.items():
            statistics.append(statistic.lower())
            for field1, field2 in fields:
                if statistic.lower() == "pk":
                    for f in [field1, field2]:
                        if f not in self.supported_fields_3D:
                            raise ValueError(f"Field {f} for statistic {statistic} not supported.")
                    fields_3D += [field1, field2]
                elif statistic.lower() == "cl" or statistic.lower() == "xi" or statistic.lower() == "cosebis":
                    for f in [field1, field2]:
                        if f not in self.supported_fields_2D:
                            raise ValueError(f"Field {f} for statistic {statistic} not supported.")
                    fields_2D += [field1, field2]
                else:
                    raise ValueError(f"Statistic {statistic} not supported.")

        # Fill in requirements if not already included
        if "xi" in statistics:
            statistics += ["pk", "cl"]
        if "cosebis" in statistics:
            statistics += ["pk", "cl"]
        if "cl" in statistics:
            statistics += ["pk"]
        if "shear" in fields_2D:
            fields_3D += ["matter"]
        if "galaxy" in fields_2D:
            fields_3D += ["matter"]
        if "y" in fields_2D:
            fields_3D += ["pressure"]

        self.statistics = list(collections.OrderedDict().fromkeys(statistics))
        self.fields_3D = list(collections.OrderedDict().fromkeys(fields_3D))
        self.fields_2D = list(collections.OrderedDict().fromkeys(fields_2D))

        self.tomographic_bin_names = tomographic_bin_names
        self.nofz_files = nofz_files

        self.parameters = parameters

        self.pipeline_has_run = False
        paths = {"HMX_PATH" : "%(KCAP_PATH)s/HMx/"}
        paths.update(module_paths)

        self.pipeline = CosmoSISPipeline(paths=paths)

        self.modules = collections.OrderedDict()

        self.add_module(("consistency", ConsistencyModule()))

        if "pk" in self.statistics:
            self.add_module(("camb", CAMBModule(transfer_function=transfer_function.lower() == "camb", 
                                                **module_configs.get("camb", {}))))
            if nonlinear_Pk.lower() == "hmx":
                self.add_module(("HMx", HMxModule(hm_mode="hmx", fields=self.fields_3D, 
                                                  transfer_function=transfer_function,
                                                  **module_configs.get("HMx", {}))))
            elif nonlinear_Pk.lower() == "hmcode":
                self.add_module(("HMx", HMxModule(hm_mode="hmcode", 
                                                  transfer_function=transfer_function,
                                                  **module_configs.get("HMx", {}))))
            else:
                raise ValueError(f"Nonlinear P(k) method {nonlinear_Pk} not supported.")

        if "cl" in self.statistics:
            self.add_module(("projection", ProjectionModule(probes=probes["cl"],
                                                            **module_configs.get("projection", {}))))
            for field in self.fields_2D:
                if field in nofz_files:
                    self.add_module((f"load_nofz_{field}", 
                                     LoadNofzModule(nofz_file=nofz_files[field],
                                                    output_section=field,
                                                    **module_configs.get("projection", {}))), before="projection")

        if "xi" in self.statistics:
            self.add_module(("cl2xi", Cl2xiModule(probes=probes["xi"],
                                                  **module_configs.get("cl2xi", {}))))

        if "cosebis" in self.statistics:
            self.add_module(("cosebis", COSEBISModule(probes=probes["cosebis"],
                                                      **module_configs.get("cosebis", {}))))

    def add_module(self, module, before=None, after=None):
        if before is None and after is None:
            self.modules[module[0]] = module[1]
        elif before is not None:
            idx = list(self.modules.keys()).index(before)
            self.modules = dict_insert(self.modules, idx, module)
        elif after is not None:
            idx = list(self.modules.keys()).index(after)
            self.modules = dict_insert(self.modules, idx, module, after=True)

    def get_pipeline_config(self):
        self.pipeline.modules = list(self.modules.values())
        ini, _ = self.pipeline.assemble_ini()
        return config_to_string(ini)

    def run_pipeline(self, parameters={}):
        self.pipeline.modules = list(self.modules.values())
        result = self.pipeline.run(parameters)
        self.pipeline_has_run = True
        return result
                
    def power_spectrum(self, probes=("matter", "matter"), k=None, z=None, kind="nonlinear"):
        if "pk" not in self.statistics:
            raise ValueError(f"TwoPoint not initialized for 3D power spectra.")
        if probes[0] not in self.fields_3D or probes[1] not in self.fields_3D:
            raise ValueError(f"TwoPoint not initialized for probes {probes}.")
        if not self.pipeline_has_run:
            self.data = self.run_pipeline(self.parameters)

        section = f"{probes[0]}_{probes[1]}_power_spectrum"
        if probes == ("matter", "matter") and not self.data.has_section(section):
            if kind.lower() == "nonlinear" and self.data.has_section("matter_power_nl"):
                section = "matter_power_nl"
            elif kind.lower() == "linear" and self.data.has_section("matter_power_lin"):
                section = "matter_power_lin"
        
        k_grid = self.data[section, "k_h"]
        z_grid = self.data[section, "z"]
        pk_grid = self.data[section, "p_k"]

        return pk_grid, k_grid, z_grid

    def extract_tomographic_bins(self, section, probe):
        twopoint = []
        bin_names = []
        for _, key in self.data.keys(section):
            if key[:3].lower() == "bin":
                bin_1, bin_2 = key.split("_")[1:3]
                twopoint.append(self.data[section, key])
                if probe[0] in self.tomographic_bin_names:
                    bin_1 = self.tomographic_bin_names[probe[0]][int(bin_1)-1]
                if probe[1] in self.tomographic_bin_names:
                    bin_2 = self.tomographic_bin_names[probe[1]][int(bin_2)-1]
                bin_names.append((bin_1, bin_2))
        return twopoint, bin_names

    def angular_power_spectrum(self, probes, ell=None,
                               return_bin_names=False):
        if "cl" not in self.statistics:
            raise ValueError(f"TwoPoint not initialized for angular power spectra.")
        if probes[0] not in self.fields_2D or probes[1] not in self.fields_2D:
            raise ValueError(f"TwoPoint not initialized for probes {probes}.")
        if not self.pipeline_has_run:
            self.data = self.run_pipeline(self.parameters)

        section = f"{probes[0]}_{probes[1]}_cl"
        if probes == ("shear", "shear") and not self.data.has_section(section):
            section = "shear_cl"
        
        cl, bin_names = self.extract_tomographic_bins(section, probes)

        cl = np.array(cl)
        ell = self.data[section, "ell"]
        if return_bin_names:
            return cl, ell, bin_names
        return cl, ell

    def correlation_function(self, probes, theta=None, kind="plus",
                             return_bin_names=False):
        if "xi" not in self.statistics:
            raise ValueError(f"TwoPoint not initialized for correlation function.")
        if probes[0] not in self.fields_2D or probes[1] not in self.fields_2D:
            raise ValueError(f"TwoPoint not initialized for probes {probes}.")
        if not self.pipeline_has_run:
            self.data = self.run_pipeline(self.parameters)

        section = f"{probes[0]}_{probes[1]}_xi"
        if probes == ("shear", "shear") and not self.data.has_section(section):
            section = f"shear_xi_{kind}"
        
        xi, bin_names = self.extract_tomographic_bins(section, probes)

        xi = np.array(xi)
        theta = self.data[section, "theta"]
        if return_bin_names:
            return xi, theta, bin_names
        return xi, theta

    def cosebis(self, probes, 
                      return_bin_names=False):
        if "cosebis" not in self.statistics:
            raise ValueError(f"TwoPoint not initialized for COSEBIS.")
        if probes[0] not in self.fields_2D or probes[1] not in self.fields_2D:
            raise ValueError(f"TwoPoint not initialized for probes {probes}.")
        if not self.pipeline_has_run:
            self.data = self.run_pipeline(self.parameters)

        section = self.modules["cosebis"].output_section_name
        cosebis, bin_names = self.extract_tomographic_bins(section, probes)

        modes = self.data[section, "cosebis_n"]
        cosebis = np.array(cosebis)
        if return_bin_names:
            return cosebis, modes, bin_names
        return cosebis, modes

# class Predict2point(cosmosis_utils.CosmoSISPipeline):
#     def __init__(self, probes=[], statistics=[], transfer_function="camb",
#                  csl_path=None, HMx_path=None, project_triad_path=None):
#         self.kcap_path = os.path.join(os.path.split(os.path.dirname(__file__))[0])
#         if csl_path is None:
#             self.csl_path = os.path.join(self.kcap_path, "cosmosis-standard-library")
#         else:
#             self.csl_path = csl_path
        
#         if HMx_path is None:
#             self.HMx_path = os.path.join(self.kcap_path, "HMx")
#         else:
#             self.HMx_path = HMx_path

#         self.cosmosis_config = self.setup_base_config()
#         self.setup_pipeline_config(self.cosmosis_config)
#         # self.setup_camb_config(self.cosmosis_config)
        
#         fields = set([p.replace("y", "pressure").replace("shear", "matter")
#                         for probe in probes 
#                             for p in probe.split("-")])
#         self.setup_HMx_config(self.cosmosis_config, fields=fields)
        
#         self.setup_extrapolate_power(self.cosmosis_config)
#         self.setup_projection_config(self.cosmosis_config)
#         self.setup_cl_filter_config(self.cosmosis_config)
#         self.setup_cl_binning_config(self.cosmosis_config)

#         # ini = 
#         # pipeline =
#         # Fast-slow split




    

#         cosmosis_config["consistency"] = {"file"              : r"%(CSL_PATH)s/utility/consistency/consistency_interface.py",}

#     def setup_camb_config(self, cosmosis_config, n_z, z_max):
#         cosmosis_config["pipeline"]["modules"] += "camb "
#         cosmosis_config["camb"] =        {"file"              : r"%(CSL_PATH)s/boltzmann/camb/camb.so",
#                                           "mode"              : "all" if self.transfer_function.lower() == "camb" else "background",
#                                           "lmax"              : "2500",
#                                           "feedback"          : "0",
#                                           "zmax"              : str(float(z_max)),
#                                           "nz"                : str(n_z),
#                                           "background_zmax"   : str(float(z_max)),
#                                           "background_nz"     : str(int(z_max*100)),}

#     def setup_HMx_config(self, cosmosis_config, n_k, k_min, k_max, n_z, z_min, z_max,
#                          fields=["dmonly"],
#                          verbose=0, hm_mode="hmcode", one_parameter_hmcode=False,
#                          ihm=None):
#         cosmosis_config["pipeline"]["modules"] += "HMx "
#         cosmosis_config["HMx"] =         {"file"              : r"%(HMX_PATH)s/lib/cosmosis_interface.so",
#                                           "fields"            : " ".join(fields),
#                                           "nk"                : str(n_k),
#                                           "kmin"              : str(float(k_min)),
#                                           "kmax"              : str(float(k_max)),
#                                           "nz"                : str(n_z),
#                                           "zmin"              : str(float(z_min)),
#                                           "zmax"              : str(float(z_max)),
#                                           "verbose"           : str(verbose),
#                                           "dimensionless_power_spectrum" : 0,
#                                           "p_lin_source"      : "external" if self.transfer_function.lower() == "camb" else "eh",
#                                           "hm_mode"           : hm_mode,
#                                           "one_parameter_hmcode" : "T" if one_parameter_hmcode else "F"}
#         if ihm is not None:
#             cosmosis_config["HMx"]["ihm"] = str(ihm)                                

#     def setup_nofz_config(self, cosmosis_config, nofz_file, output_section="nz_shear"):
#         if isinstance(nofz_file, str):
#             filepath = nofz_file
#         elif isinstance(nofz_file, (list, tuple)):
#             filepath = " ".join(nofz_file)
#         else:
#             raise ValueError("nofz_file needs to be str or list of str.")
#         cosmosis_config["pipeline"]["modules"] += "load_nz "
#         cosmosis_config["load_nz"] =     {"file"              : r"%(CSL_PATH)s/number_density/load_nz/load_nz.py",
#                                           "filepath"          : filepath,
#                                           "histogram"         : "True",
#                                           "output_section"    : output_section,}

#     def setup_projection_config(self, cosmosis_config, probes,
#                                 n_ell=400, ell_min=0.1, ell_max=5e4, verbose=True):
#         need_tSZ = False
#         for probe in probes:
#             if "y" in probe:
#                 need_tSZ = True
#         if need_tSZ:
#             module_file = r"%(PROJECT_TRIAD_PATH)s/projection_tSZ/project_2d_tSZ.py"
#         else:
#             module_file = r"%(CSL_PATH)s/cosmosis-standard-library/structure/projection/project_2d.py"
#         cosmosis_config["pipeline"]["modules"] += "projection "
#         cosmosis_config["projection"] =  {"file"              : module_file,
#                                           "ell_min"           : str(float(ell_min)),
#                                           "ell_max"           : str(float(ell_max)),
#                                           "n_ell"             : str(n_ell),
#                                           "verbose"           : "T" if verbose else "F",
#                                           "get_kernel_peaks"  : "F"}
#         for probe in probes:
#             cosmosis_config["projection"][probe] = probe.replace("cmbkappa", "k")

#     def setup_cl_filter_config(self, cosmosis_config, probes, y_filter_fwhm):
#         cosmosis_config["pipeline"]["modules"] += "filter_Cls "
#         sections = []
#         powers = []
#         for probe in probes:
#             p1, p2 = probe.split("-")
#             if p1 == "y" or p2 == "y":
#                 sections.append(f"{probe.replace('-', '_')}_cl")
#                 if p1 == "y" and p2 == "y":
#                     powers.append("2")
#                 else:
#                     powers.append("1")

#         cosmosis_config["filter_Cls"] =           {"file"              : r"%(PROJECT_TRIAD_PATH)s/tools/filter_cls.py",
#                                                    "filter"            : "gaussian",
#                                                    "fwhm"              : str(y_filter_fwhm),
#                                                    "sections"          : " ".join(sections),
#                                                    "powers"            : " ".join(powers),}

#     def setup_cl_binning_config(self, cosmosis_config, n_bin, ell_min, ell_max, logspaced=True):
#         cosmosis_config["pipeline"]["modules"] += "bin_Cls "
#         cosmosis_config["bin_Cls"] =              {"file"              : r"%(PROJECT_TRIAD_PATH)s/tools/bin_cls.py",
#                                                    "ell_min"           : str(ell_min),
#                                                    "ell_max"           : str(ell_max),
#                                                    "n_bin"             : str(n_bin),
#                                                    "logspaced"         : "T" if logspaced else "F",}
        
