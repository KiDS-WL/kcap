from cosmosis.datablock import option_section, names

# Modified from add_instrinsic
def setup(options):
    do_position_position = options.get_bool(option_section, "position-position", False)
    do_position_shear = options.get_bool(option_section, "position-shear", True)

    if do_position_position:
        raise NotImplementedError("add_magnification only does GGL for now.")

    suffix = options.get_string(option_section, "suffix", "")

    if suffix:
        suffix = "_" + suffix

    sec_names = {
        "galaxy_shear": "galaxy_shear_cl" + suffix,
        "magnification_shear": "magnification_shear_cl" + suffix,
    }   

    return do_position_position, do_position_shear, sec_names


def execute(block, config):
    do_position_position, do_position_shear, sec_names = config

    galaxy_shear = sec_names['galaxy_shear']
    magnification_shear = sec_names['magnification_shear']
    
    if do_position_shear:
        nbin_shear = block[galaxy_shear, 'nbin_b']
        nbin_pos = block[galaxy_shear, 'nbin_a']

    if do_position_shear:
        for i in range(nbin_pos):
            for j in range(nbin_shear):
                bin_ij = f"bin_{i+1}_{j+1}"
                block[galaxy_shear, bin_ij] += block[magnification_shear, bin_ij]



    return 0
