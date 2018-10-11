from __future__ import print_function
from builtins import range
import numpy as np
from cosmosis.datablock import option_section, names as section_names
from cosmosis.datablock.cosmosis_py.errors import BlockWrongValueType

def setup(options):
    # only one parameter - filepath
    filename = options[option_section, "filepath"]
    des_fmt = options.get_bool(option_section, "des_fmt", default=False)
    histogram = options.get_bool(option_section, "histogram", default=False)

    output_section = options.get_string(
        option_section, "output_section", default=section_names.wl_number_density)
    single_bin = options.get_int(option_section, "single_bin", default=-666)
    upsampling = options.get_int(option_section, "upsampling", default=1)

    try:
        bins = options[option_section, "bins"]
        if isinstance(bins, int):
            single_bin = bins
    except:
        bins = None

    if " " in filename:
        # Multiple files.
        filenames = [f.strip() for f in filename.split(" ")]
        data = []
        for f in filenames:
            data.append(np.loadtxt(f).T)
        # Check that z values are consistent
        z = data[0][0]
        if not all([np.allclose(z, d[0]) for d in data]):
            raise ValueError("Mismatch of z column in n(z) files.")
        data_full = np.concatenate([z[None,:]]+[d[1:] for d in data], axis=0)
    else:
        data_full = np.loadtxt(filename).T

    if des_fmt:
        z = 0.5 * (data_full[0] + data_full[1])
        nz = len(z)
        nbin = len(data_full) - 3
        n_of_z = data_full[2:-1]
    else:
        nz = len(data_full[0])
        z = data_full[0]
        if single_bin != -666:
            n_of_z = data_full[single_bin]
            nbin = 1
        else:
            if bins is None:
                nbin = len(data_full) - 1
                n_of_z = data_full[1:]
            else:
                nbin = len(bins)
                n_of_z = data_full[bins]

    if histogram:
        # in this case the sample z values are lower edges of
        # histogram bins.  So to turn them into samples we need to
        # shift everything.  This assumes equal sized bins
        dz = (z[1] - z[0]) / 2.0
        print("n(z) set to histogram mode. Bin centers are %f higher than edges." % dz)
        z += dz

    # check first z is zero, if not add some
    if z[0] > 0.00000001:
        z_new = np.zeros(len(z) + 1)
        z_new[1:] = z
        n_of_z_new = np.zeros((nbin, len(z) + 1))
        n_of_z_new[:, 1:] = n_of_z
        z, n_of_z = z_new, n_of_z_new

    if upsampling > 1:
        print("Upsampling z by factor {}".format(upsampling))
        z_new = np.linspace(0.0, z[-1], len(z) * upsampling)
        sample_bin = np.digitize(z_new, z) - 1
        n_of_z_new = np.zeros((nbin, len(z_new)))
        for i in range(nbin):
            n_of_z_new[i][:] = n_of_z[i][sample_bin]
        z, n_of_z = z_new, n_of_z_new

    # Normalize n(z)
    for col in n_of_z:
        norm = np.trapz(col, z)
        col /= norm

    print("Found %d samples and %d bins in redshift in file %s" % (nbin, nz, filename))
    return (nz, nbin, z, n_of_z, output_section)


def execute(block, config):
    (nz, nbin, z, n_of_z, output_section) = config

    block[output_section, 'nz'] = nz
    block[output_section, 'nbin'] = nbin
    block[output_section, 'z'] = z

    for (bin, bin_n_of_z) in enumerate(n_of_z):
        name = "bin_%d" % (bin + 1)
        block[output_section, name] = bin_n_of_z

    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
