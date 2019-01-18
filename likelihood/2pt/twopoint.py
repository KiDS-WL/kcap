from __future__ import print_function
from builtins import zip
from builtins import range
from builtins import object
from astropy.io import fits
import astropy.units
from astropy.table import Table
from enum34 import Enum
import numpy as np
import copy
# FITS header keyword indicating 2-point data
TWOPOINT_SENTINEL = "2PTDATA"
NZ_SENTINEL = "NZDATA"
COV_SENTINEL = "COVDATA"
window_types = ["SAMPLE", "CLBP", "LOG_MID"]
#LOG_MID interpolates at log mid of bin

# Please do not add things to this list
ANGULAR_UNIT_TYPES = [
    astropy.units.arcsec,
    astropy.units.arcmin,
    astropy.units.rad,
    astropy.units.deg,
]
ANGULAR_UNITS = {unit.name: unit for unit in ANGULAR_UNIT_TYPES}

def sample_cov(xi_arrays, mode='full'):
    """mode should be full, subsample or jk"""
    nsample, npoints = xi_arrays.shape
    Cov = np.zeros((npoints, npoints))
    Corr = np.zeros_like(Cov)

    xi_mean = np.mean(xi_arrays, axis=0)
    for i in range(npoints):
        for j in range(npoints):
            Cov[i, j] = np.sum((xi_arrays[:, i] - xi_mean[i])
                               * (xi_arrays[:, j] - xi_mean[j])) / nsample
    # This is the covariance in a patch of size A/nsample. Assume Cov ~ 1/A so Cov=Cov/nsample
    if mode == 'subsample':
        Cov /= nsample
    elif mode == "jk":
        Cov *= (nsample - 1)
    for i in range(npoints):
        for j in range(npoints):
            Corr[i, j] = Cov[i, j] / np.sqrt(Cov[i, i] * Cov[j, j])
    return Cov, Corr


class Types(Enum):
    """
    This is an enumeration - a list of possible values with names and code values
    that we can use in FITS headers.

    It enumerates the different quantities that a two-point measurement can correlate.
    For example, CMB T,E,B, lensing E and B, galaxy position and magnification.

    It specifies the quantity and whether we are in Fourier space or real space.

    One special case is xi_{-} and xi_{+} in galaxy shear.  These values are already
    correlations of combinations of E and B modes. We denote this as xi_{++} and xi_{--} here.

    """
    galaxy_position_fourier = "GPF"
    galaxy_shear_emode_fourier = "GEF"
    galaxy_shear_bmode_fourier = "GBF"
    galaxy_position_real = "GPR"
    galaxy_shear_plus_real = "G+R"
    galaxy_shear_minus_real = "G-R"
    cmb_kappa_real = "CKR"

    @classmethod
    def lookup(cls, value):
        for T in cls:
            if T.value == value:
                return T


def dummy_kernel(name):
    return NumberDensity(name, np.zeros(10), np.zeros(10), np.zeros(10), [np.zeros(10)])


class NumberDensity(object):
    """
    This class contains n(z) information for a collection of numbered redshift bins.
    It is expected to be used for a single sample type (e.g. galaxy sample) split into
    tomographic bins rather than the more complex cases where there are multiple different
    quantities (e.g. two galaxy samples).

    Since the main use case for this is photometric redshifts, and photo-z codes typically
    produce histogram type data sets (that is, they look like step functions between each bin),
    that is what this form assumes.
    """

    def __init__(self, name, zlow, z, zhigh, nzs, ngal=None, sigma_e=None):
        self.name = name
        self.zlow = zlow
        self.z = z
        self.zhigh = zhigh
        self.nbin = len(nzs)
        if self.nbin > 0:
            self.nsample = len(nzs[0])
        else:
            self.nsample = 0
        self.nzs = nzs
        self.ngal = ngal
        self.sigma_e = sigma_e

    @classmethod
    def from_fits(cls, extension):
        # load in the n(z) data
        data = extension.data
        header = extension.header

        z = data['Z_MID']
        zlow = data['Z_LOW']
        zhigh = data['Z_HIGH']
        i = 1
        name = 'BIN{}'.format(i)
        nzs = []
        sigma_e = []
        ngal = []
        while name in data.names:
            nz = data[name]
            nzs.append(nz)
            ngal.append(header.get('NGAL_{}'.format(i)))
            sigma_e.append(header.get('SIG_E_{}'.format(i)))
            i += 1
            name = 'BIN{}'.format(i)

        if all(x is None for x in sigma_e):
            sigma_e = None
        else:
            assert not any(x is None for x in sigma_e), "Please specify all or none of the SIG_E_ in your n(z) section {}".format(
                extension.name)
            sigma_e = np.array(sigma_e)

        if all(x is None for x in ngal):
            ngal = None
        else:
            assert not any(x is None for x in ngal), "Please specify all or none of the NGAL in your n(z) section {}".format(
                extension.name)
            ngal = np.array(ngal)

        N = cls(extension.name, zlow, z, zhigh,
                nzs, ngal=ngal, sigma_e=sigma_e)

        return N

    def to_fits(self):
        header = fits.Header()
        header[NZ_SENTINEL] = True
        header['EXTNAME'] = self.name

        columns = [
            fits.Column(name='Z_LOW', array=self.zlow, format='D'),
            fits.Column(name='Z_MID', array=self.z, format='D'),
            fits.Column(name='Z_HIGH', array=self.zhigh, format='D'),
        ]

        for i in range(self.nbin):
            name = "BIN{}".format(i + 1)
            columns.append(fits.Column(
                name=name, array=self.nzs[i], format='D'))

        if self.sigma_e is not None:
            for i in range(self.nbin):
                name = "SIG_E_{}".format(i + 1)
                header[name] = self.sigma_e[i]

        if self.ngal is not None:
            for i in range(self.nbin):
                name = "NGAL_{}".format(i + 1)
                header[name] = self.ngal[i]

        extension = fits.BinTableHDU.from_columns(columns, header=header)
        return extension

    @classmethod
    def from_block(cls, block, section_name, output_name=None):
        """Load kernel from cosmosis datablock"""
        if output_name==None:
            output_name = section_name
        if block.has_value( section_name, "z_mid" ):
            z_mid = block[section_name, "z_mid"]
            z_low = block[section_name, "z_low"]
            z_high = block[section_name, "z_high"]
        else:
            z_mid = block[section_name, "z"]
            dz = z_mid[2]-z_mid[1]
            z_low = z_mid - 0.5 * dz
            z_high = z_mid + 0.5 * dz

        nzs = []
        for i in range(1,99999):
            bin_name = "bin_%d"%i
            if block.has_value( section_name, bin_name ):
                nzs.append( block[section_name, bin_name] )
            else:
                break

        N = cls( output_name, z_low, z_mid, z_high,
                nzs )
        return N



class SpectrumMeasurement(object):
    def __init__(self, name, bins, types, kernels, windows, angular_bin, value,
                 angle=None, error=None, angle_unit=None, metadata=None,
                 angle_min=None, angle_max=None):
        """metadata is a dictionary which will get added to the fits header"""
        self.name = name
        self.bin1, self.bin2 = bins
        self.bin_pairs = self.get_bin_pairs()  # unique bin pairs
        self.type1, self.type2 = types
        self.kernel1, self.kernel2 = kernels
        self.angular_bin = angular_bin
        self.angle = angle
        self.angle_min = angle_min
        self.angle_max = angle_max
        self.value = value
        if windows in window_types:
            self.windows = windows
        else:
            raise TypeError("window type %s not recognised" % windows)
        self.error = error
        self.metadata = metadata
        if self.is_real_space():
            #angle is real
            msg = "Files with real-space units must specify units as one of: {}".format(
                list(ANGULAR_UNITS.keys()))
            assert angle_unit in ANGULAR_UNITS,  msg
        self.angle_unit = angle_unit

    def get_bin_pairs(self):
        all_pairs = list(zip(self.bin1, self.bin2))
        unique_pairs = []
        for p in all_pairs:
            if p not in unique_pairs:
                unique_pairs.append(p)
        return unique_pairs

    def is_real_space(self):
        return self.type1.value.endswith("R") or self.type2.value.endswith("R")

    def convert_angular_units(self, unit):
        if not self.is_real_space():
            raise ValueError(
                "Two point spectrum has no units to convert; it is in Fourier space")

        if self.windows not in ["SAMPLE", "CLBP"]:
            raise NotImplementedError(
                "Need to write code to transform window function units")
        old_unit = ANGULAR_UNITS[self.angle_unit]
        new_unit = ANGULAR_UNITS[unit]

        print("Converting angle units of {} from {} -> {} (factor {})".format(
            self.name, old_unit, new_unit, old_unit.to(new_unit)))
        angle_with_units = self.angle * old_unit
        self.angle = angle_with_units.to(new_unit).value
        self.angle_unit = unit
        if self.angle_min is not None:
            self.angle_min = ( self.angle_min * old_unit).to(new_unit).value
            self.angle_max = ( self.angle_max * old_unit).to(new_unit).value

    def apply_mask(self, mask):
        """mask is a boolean array which True for elements to be kept"""
        self.bin1 = self.bin1[mask]
        self.bin2 = self.bin2[mask]
        self.angular_bin = self.angular_bin[mask]
        self.angle = self.angle[mask]
        if self.angle_min is not None:
            self.angle_min = self.angle_min[mask]
            self.angle_max = self.angle_max[mask]
        self.value = self.value[mask]
        if self.error is not None:
            self.error = self.error[mask]

    def cut_bin_pair(self, bin_pair, complain=False):
        """Cut a full bin pair. If complain is True,
        raise a ValueError if the bin pair is not found"""
        b1, b2 = bin_pair
        mask = (self.bin1==b1)*(self.bin2==b2)
        if (mask.sum()==0) and complain:
            raise ValueError("""You tried to cut bin pair %d,%d from spectrum named %s 
                             but it doesn't exist"""%(b1,b2,self.name))
        self.apply_mask(~mask)

    def auto_bins(self):
        return self.bin1 == self.bin2

    def __len__(self):
        return len(self.value)

    def nbin(self):
        return np.max([self.bin1.max(), self.bin2.max()])

    def get_pair(self, bin1, bin2):
        w = (self.bin1 == bin1) & (self.bin2 == bin2)
        return self.angle[w], self.value[w]

    def get_pair_mask(self, bin1, bin2):
        w = (self.bin1 == bin1) & (self.bin2 == bin2)
        return w

    def get_error(self, bin1, bin2):
        if self.error is None:
            return None
        w = (self.bin1 == bin1) & (self.bin2 == bin2)
        return self.error[w]

    @classmethod
    def from_fits(cls, extension, covmat_info=None):
        name = extension.name
        # determine the type of the quantity involved
        type1 = Types.lookup(extension.header['QUANT1'])
        type2 = Types.lookup(extension.header['QUANT2'])

        # and the name of the kernels and the window function
        # extensions
        kernel1 = extension.header['KERNEL_1']
        kernel2 = extension.header['KERNEL_2']
        windows = extension.header['WINDOWS']

        if windows not in window_types:
            raise NotImplementedError(
                "Have not yet coded window functions for angular bins")

        # Now load the data
        data = extension.data
        bin1 = data['BIN1']
        bin2 = data['BIN2']
        angular_bin = data['ANGBIN']
        value = data['VALUE']
        if "ANG" in data.names:
            angle = data['ANG']
            ang_index = data.names.index("ANG")
            angle_unit = extension.header.get('TUNIT{}'.format(ang_index + 1))
            if angle_unit is not None:
                angle_unit = angle_unit.strip()
        else:
            angle = None
            angle_unit = None

        # Load a chunk of the covariance matrix too if present.
        if covmat_info is None:
            error = None
        else:
            error = covmat_info.get_error(name)

        return SpectrumMeasurement(name, (bin1, bin2), (type1, type2), (kernel1, kernel2), windows,
                                   angular_bin, value, angle, error, angle_unit=angle_unit)

    def to_fits(self):
        header = fits.Header()
        header[TWOPOINT_SENTINEL] = True
        header['EXTNAME'] = self.name
        header['QUANT1'] = self.type1.value
        header['QUANT2'] = self.type2.value
        header['KERNEL_1'] = self.kernel1
        header['KERNEL_2'] = self.kernel2
        header['WINDOWS'] = self.windows  # NOT YET CODED ANYTHING ELSE
        header['N_ZBIN_1'] = len(np.unique(self.bin1))
        header['N_ZBIN_2'] = len(np.unique(self.bin2))
        if self.metadata is not None:
            # Check metadata doesn't share any keys with the stuff that's already in the header
            assert set(self.metadata.keys()).isdisjoint(list(header.keys()))
            for key, val in list(self.metadata.items()):
                header[key] = val
        header['N_ANG'] = len(np.unique(self.angular_bin))

        columns = [
            fits.Column(name='BIN1', array=self.bin1, format='K'),
            fits.Column(name='BIN2', array=self.bin2, format='K'),
            fits.Column(name='ANGBIN', array=self.angular_bin, format='K'),
            fits.Column(name='VALUE', array=self.value, format='D'),
        ]
        if self.angle is not None:
            if self.windows == "SAMPLE":
                columns.append(fits.Column(
                    name='ANG', array=self.angle, format='D', unit=self.angle_unit))
            if self.windows == "CLBP":
                columns.append(fits.Column(
                    name='ANG', array=self.angle, format='2K', unit=self.angle_unit))

        extension = fits.BinTableHDU.from_columns(columns, header=header)
        return extension

class CovarianceMatrixInfo(object):
    """Encapsulate a covariance matrix and indices in it"""

    def __init__(self, name, names, lengths, covmat):
        super(CovarianceMatrixInfo, self).__init__()
        self.name = name
        self.names = names
        self.lengths = lengths
        self.starts = [0]
        for i, l in enumerate(self.lengths[:-1]):
            self.starts.append(l + self.starts[i])
        self.covmat = covmat
        self.diagonal = covmat.diagonal()

    def get_error(self, name):
        i = self.names.index(name)
        start = self.starts[i]
        end = start + self.lengths[i]
        return self.diagonal[start:end]**0.5

    def to_fits(self):
        header = fits.Header()
        header[COV_SENTINEL] = True
        header['EXTNAME'] = self.name
        for i, (start_index, name) in enumerate(zip(self.starts, self.names)):
            header['STRT_{}'.format(i)] = start_index
            header['NAME_{}'.format(i)] = name
        extension = fits.ImageHDU(data=self.covmat, header=header)
        return extension

    @classmethod
    def from_fits(cls, extension):
        cov_name = extension.name
        covmat = extension.data
        header = extension.header
        i = 0
        measurement_names = []
        start_indices = []
        while True:
            name_card = 'NAME_{}'.format(i)
            if name_card not in header:
                break
            measurement_names.append(header[name_card])
            start_index_card = 'STRT_{}'.format(i)
            start_indices.append(header[start_index_card])
            i += 1
        lengths = []
        current_length = 0
        # this only works if more than one spectrum
        if len(start_indices) > 1:
            for start, end in zip(start_indices[:-1], start_indices[1:]):
                lengths.append(end - start)
            if start_indices:
                lengths.append(covmat.shape[0] - start_indices[-1])
        else:
            lengths.append(covmat.shape[0])
        return cls(cov_name, measurement_names, lengths, covmat)

    @classmethod
    def from_spec_lists(cls, spec_lists, cov_name, mode='full'):
        """Often the covariance will be computed by measuring the statistic(s) in question
        on many simulated realisations of the dataset. This function takes a list of such 
        measurements, *each one a list SpectrumMeasurement objects*, computes the mean and covariance, and
        returns the mean as a list of SpectrumMeasurements, and the covariance as a CovarianceMatrixInfo
        object. mode should be one of full, subsample or jackknife"""
        # first check that spec_lists is a list of lists, and that there are at least 2
        print('spec_lists', spec_lists)
        try:
            spec_lists[1][0]
        except Exception as e:
            print("spec_lists should be a list of lists with at least two elements")
            raise(e)

        # Get measurement names and lengths from first element of spec_lists
        num_spec = len(spec_lists[0])
        names = [s.name for s in spec_lists[0]]
        lengths = [len(s.value) for s in spec_lists[0]]

        # Now loop through all realisations, building up list of numpy arrays of raw measurements
        n_real = len(spec_lists)
        spec_arrays = []
        for i_real in range(n_real):
            spec_array = []
            # check this realisation has right number,type,length of spectra
            for i_spec in range(num_spec):
                assert spec_lists[i_real][i_spec].name == names[i_spec]
                assert len(spec_lists[i_real][i_spec].value) == lengths[i_spec]
                spec_array += list(spec_lists[i_real][i_spec].value)
            spec_arrays.append(np.array(spec_array))

        # Now compute covariance
        spec_arrays = np.array(spec_arrays)
        cov_values, _ = sample_cov(spec_arrays, mode=mode)
        mean_spec_values = np.mean(spec_arrays, axis=0)

        # make list of mean 2pt specs by copying spec_lists[0] and replacing value column
        mean_spec = copy.copy(spec_lists[0])
        index_start = 0
        for i_spec in range(num_spec):
            end = index_start + lengths[i_spec]
            inds = np.arange(index_start, index_start + lengths[i_spec])
            index_start = end
            mean_spec[i_spec].value = mean_spec_values[inds]

        return cls(cov_name, names, lengths, cov_values), mean_spec


class TwoPointFile(object):
    def __init__(self, spectra, kernels, windows, covmat_info):
        if windows is None:
            windows = {}
        self.spectra = spectra
        self.kernels = kernels
        self.windows = windows
        self.covmat_info = covmat_info
        if covmat_info:
            #self.covmat = covmat_info.covmat
            self.covmat = self.get_cov_start()
        else:
            self.covmat = None

    def get_spectrum(self, name):
        spectra = [spectrum for spectrum in self.spectra if spectrum.name == name]
        n = len(spectra)
        if n == 0:
            raise ValueError("Spectrum with name %s not found in file" % name)
        elif n > 1:
            raise ValueError(
                "Multiple spectra with name %s found in file" % name)
        else:
            return spectra[0]

    def get_kernel(self, name):
        kernels = [kernel for kernel in self.kernels if kernel.name == name]
        n = len(kernels)
        if n == 0:
            raise ValueError("Kernel with name %s not found in file" % name)
        elif n > 1:
            raise ValueError(
                "Multiple kernel with name %s found in file" % name)
        else:
            return kernels[0]

    def _mask_covmat(self, masks):
        # Also cut down the covariance matrix accordingly
        if self.covmat is not None:
            mask = np.concatenate(masks)
            self.covmat = self.covmat[mask, :][:, mask]

    def mask_bad(self, bad_value):
        "Go through all the spectra masking out data where they are equal to bad_value"
        masks = []
        # go through the spectra and covmat, masking out the bad values.
        for spectrum in self.spectra:
            # nb this will not work for NaN!
            mask = (spectrum.value != bad_value)
            spectrum.apply_mask(mask)
            print("Masking {} values in {}".format(mask.size - mask.sum(), spectrum.name))
            # record the mask vector as we will need it to mask the covmat
            masks.append(mask)
        if masks:
            self._mask_covmat(masks)

    def mask_cross(self):
        masks = []
        for spectrum in self.spectra:
            mask = spectrum.auto_bins()
            spectrum.apply_mask(mask)
            print("Masking {} cross-values in {}".format(mask.size - mask.sum(), spectrum.name))
            masks.append(mask)
        if masks:
            self._mask_covmat(masks)

    def mask_scales(self, cuts={}, bin_cuts=[]):
        masks = []
        print()
        for spectrum in self.spectra:
            mask = np.ones(len(spectrum), dtype=bool)
            for b1, b2 in spectrum.bin_pairs:
                w_full = np.where((spectrum.bin1 == b1) &
                                  (spectrum.bin2 == b2))[0]
                if (spectrum.name, b1, b2) in bin_cuts:
                    print("Removing {} bin ({},{}) altogether.".format(spectrum.name, b1, b2))
                    mask[w_full] = False
                    continue

                cut = cuts.get((spectrum.name, b1, b2))
                if cut is None:
                    print("No cut specified for {} bin ({},{})".format(spectrum.name, b1, b2))
                    continue

                # Actually do the cut
                ang_min, ang_max = cut
                w = np.where((spectrum.bin1 == b1) & (spectrum.bin2 == b2) &
                             ((spectrum.angle < ang_min) | (spectrum.angle > ang_max)))[0]

                print("Cutting {} bin pair ({},{}) to angle range ({} - {}) : this removes {} values out of {}".format(
                    spectrum.name, b1, b2, ang_min, ang_max, len(w), len(w_full)))

                mask[w] = False
            masks.append(mask)
            spectrum.apply_mask(mask)
            print()

        if masks:
            self._mask_covmat(masks)

    def mask_scale(self, spectra_to_cut, min_scale=-np.inf, max_scale=np.inf):
        masks = []
        # go through the spectra and covmat, masking out the bad values.
        for spectrum in self.spectra:
            mask = np.ones(len(spectrum.value), dtype=bool)
            if (spectra_to_cut != "all") and (spectrum.name not in spectra_to_cut):
                masks.append(mask)
            else:
                # nb this will not work for NaN!
                mask = (spectrum.angle > min_scale) & (
                    spectrum.angle < max_scale)
                spectrum.apply_mask(mask)
                print("Masking {} values in {} because they had ell or theta outside ({},{})".format(mask.size - mask.sum(), spectrum.name, min_scale, max_scale))
                # record the mask vector as we will need it to mask the covmat
                masks.append(mask)

        if masks:
            self._mask_covmat(masks)

    def choose_data_sets(self, data_sets):
        """Strip out any data sets not in the given list."""
        data_sets = [d.lower() for d in data_sets]
        mask = []
        use = []
        for spectrum in self.spectra:
            if spectrum.name.lower() in data_sets:
                use.append(True)
                mask.append(np.ones(spectrum.bin1.size, dtype=bool))
            else:
                use.append(False)
                mask.append(np.zeros(spectrum.bin1.size, dtype=bool))
        for data_set in data_sets:
            if not any(spectrum.name.lower() == data_set for spectrum in self.spectra):
                raise ValueError(
                    "Data set called {} not found in two-point data file.".format(data_set))
        self.spectra = [s for (u, s) in zip(use, self.spectra) if u]
        if self.covmat is not None:
            mask = np.concatenate(mask)
            self.covmat = self.covmat[mask, :][:, mask]

    def get_cov_start(self):

        # This gets the covariance array in the right order (before any scale cuts etc.)
        cov = self.covmat_info.covmat
        if self.covmat_info.names == [spec.name for spec in self.spectra]:
            # Ordering is already ok
            return cov
        print("Covariance matrix is not in the same order as the 2pt measurement extensions...doing some damn fiddly")
        print("re-ordering, if I screw it up it's your fault for not putting your covariance in the right order")
        cov_starts = self.covmat_info.starts
        cov_lengths = self.covmat_info.lengths
        cov_names = self.covmat_info.names
        cov_ends = [cov_lengths[0]]
        for i in range(len(cov_names) - 1):
            cov_ends.append(cov_ends[i] + cov_lengths[i + 1])
        # print 'cov_lengths',cov_lengths
        # print 'cov_starts',cov_starts
        # print 'cov_ends',cov_ends
        assert cov_ends[-1] == cov.shape[0]

        total_l = 0
        spec_inds = []
        spec_names = [spec.name for spec in self.spectra]
        for spectrum in self.spectra:
            spec_inds.append(cov_names.index(spectrum.name))
            total_l += cov_lengths[cov_names.index(spectrum.name)]
        cov_out = np.zeros((total_l, total_l))
        start_i = 0
        # print spec_names
        # print spec_inds

        for ti, ind_i in zip(spec_names, spec_inds):
            start_j = 0
            for tj, ind_j in zip(spec_names, spec_inds):
                cov_out[start_i:start_i + cov_lengths[ind_i], start_j:start_j + cov_lengths[ind_j]
                        ] = cov[cov_starts[ind_i]:cov_ends[ind_i], cov_starts[ind_j]:cov_ends[ind_j]]
                start_j += cov_lengths[ind_j]
            start_i += cov_lengths[ind_i]
        return cov_out

    def to_fits(self, filename, overwrite=False):
        hdus = [fits.PrimaryHDU()]

        if self.covmat_info is not None:
            hdus.append(self.covmat_info.to_fits())

        for spectrum in self.spectra:
            if spectrum.windows not in window_types:
                raise NotImplementedError(
                    "Sorry - not yet coded general case with ell/theta window functions")
            hdus.append(spectrum.to_fits())

        if self.kernels is not None:
            for kernel in self.kernels:
                hdus.append(kernel.to_fits())

        hdulist = fits.HDUList(hdus)
        hdulist.writeto(filename, overwrite=overwrite)

    @classmethod
    def from_fits(cls, filename, covmat_name="COVMAT"):
        fitsfile = fits.open(filename)
        spectra = []
        kernels = []
        windows = {}

        # Load the covariance matrix from the file, typically stored as
        # image data.
        # It turns out it is more conventient to load this first.
        # We can also use no covmat at all.
        if covmat_name is None:
            covmat_info = None
            covmat = None
        else:
            extension = fitsfile[covmat_name]
            covmat_info = CovarianceMatrixInfo.from_fits(extension)

        # First load all the spectra in the file
        # Spectra are indicated by the "2PTDATA" keyword
        # in the header being "T"
        for extension in fitsfile:
            if extension.header.get(TWOPOINT_SENTINEL):
                spectra.append(SpectrumMeasurement.from_fits(
                    extension, covmat_info))

        # Each spectrum needs kernels, usually some n(z).
        # These were read from headers when we loaded the 2pt data above above.
        # Now we are loading those kernels into a dictionary
        for extension in fitsfile:
            if extension.header.get(NZ_SENTINEL):
                kernels.append(NumberDensity.from_fits(extension))

        # We might also require window functions W(ell) or W(theta). It's also possible
        # that we just use a single sample value of ell or theta instead.
        # If the spectra required it (according to the header) then we also
        # need to load those in.
        for spectrum in spectra:
            if spectrum.windows not in window_types and spectrum.windows not in windows:
                windows[spectrum.windows] = cls._windows_from_fits(
                    fitsfile[windows])

        # return a new TwoPointFile object with the data we have just loaded
        return cls(spectra, kernels, windows, covmat_info)

    @classmethod
    def _windows_from_fits(cls, extension):
        raise NotImplementedError("non-sample window functions in ell/theta")
