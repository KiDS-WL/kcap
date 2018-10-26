import os
import glob

import numpy as np

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