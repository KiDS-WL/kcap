"""
This is a template for an example script showing how you 
could add a covariance matrix to an existing file.
"""
import twopoint
import numpy as np
import sys


twopoint_filename = sys.argv[1]
covmat_filename = sys.argv[2]
new_filename = sys.argv[3]

print("Loading 2pt data in {}".format(twopoint_filename))
data = twopoint.TwoPointFile.from_fits(twopoint_filename, covmat_name=None)

print("Loading text covariance file from {}".format(covmat_filename))
covmat = np.loadtxt(covmat_filename)

print("Now you need to make sure that the ordering of the covariance matrix")
print("matches the ordering in the text file!")
print("Expecting this order:")
print("#Name Bin1 Bin2 Angle")
for s in data.spectra:
	for b1, b2, ang in zip(s.bin1, s.bin2, s.angle):
		print(s.name, b1, b2, ang)


replace_this_with_some_reordering_if_needed()

# The ordering of the 
names = [s.name for s in data.spectra]
lengths = [len(s) for s in data.spectra]
n = sum(lengths)

assert covmat.shape==(n,n)

data.covmat_info = twopoint.CovarianceMatrixInfo("COVMAT", names, lengths, covmat)
data.to_fits(new_filename)
