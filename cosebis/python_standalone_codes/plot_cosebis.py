from   array import array
import matplotlib.pylab as plt
import matplotlib.pyplot
from   matplotlib.font_manager import FontProperties
from   matplotlib.ticker import ScalarFormatter
import numpy as np
import math, os
import glob, pickle
import matplotlib.image as mpimg
from numpy.linalg import inv
from numpy import linalg as LA
from matplotlib.patches import Rectangle
from scipy import interpolate
from cutCov import cutCov



arcmin=np.pi/180./60.

nMaximum=20
nMin=0
nMax=5
tmin=0.5
tmax=300.0

nBins=5
nPairs= int((nBins+1)*nBins/2)

thetamin='%.2f' % tmin
thetamax='%.2f' % tmax
thetaRange=thetamin+'-'+thetamax  
thetamax='%.0f' % tmax
thetamin='%.1f' % tmin
thetaTitle=r'$['+thetamin+"', "+thetamax+"']$"
thetaFileName=thetamin+'_'+thetamax

FolderName="../../cosmosis_in_out/KV450_analysis/"
FolderName_out=FolderName+"outputs/cosebis_HHi/output_Hildebrandt_lnAs_test_withScaleCuts_nMax"+str(nMax)+"/"
output_folder ="plots/"

filename=FolderName+"inputs/Covariance_nMaximum_20_0.50_300.00_nBins5_sigma_m_0.0200.ascii"
file=open(filename)
cov_in=np.loadtxt(file,comments='#')


filename=FolderName+"inputs/Cov_SSC_thKV450_nell_500_20_0.50-300.00_nBins_5.ascii"
file=open(filename)
cov=np.loadtxt(file,comments='#')
# cut the range needed
cov_in+=cov

filename=FolderName+"inputs/Covariance_nMaximum_20_0.50_300.00_nBins5_NoiseJustForNoise.ascii"
file=open(filename)
cov=np.loadtxt(file,comments='#')
# cut the range needed
cov_in+=cov

cov_in=cutCov(cov_in,nPairs,nMaximum,nMin,nMax)


filename=FolderName_out+"/likelihood_en/covariance.txt"
file=open(filename)
CovEn_scalecut=np.loadtxt(file,comments='#')

FolderName_out_alreadycut=FolderName+"outputs/cosebis_HHi/output_Hildebrandt_lnAs_test_withScaleCuts/"
filename=FolderName_out+"/likelihood_en/covariance.txt"
file=open(filename)
CovEn_alreadycut=np.loadtxt(file,comments='#')

if((abs(CovEn_scalecut/cov_in-1)<1e-10).all()):
	print("in and out covs are the same")
else:
	print("covariances are not the same!")

exit(1)

for b in range(nBands):
	iband=b+1
	ellRange='%.2f' % l_min_vec[b]+"-"+ '%.2f' % l_max_vec[b]
	Wfilename=FolderName+"W_noAp_"+str(bessel_order)+"_"+ellRange+"_"+thetaRange+".table"	
	file=open(Wfilename)
	W=np.loadtxt(file,comments='#')
# 
	Wfilename=FolderName+"W_ap0.50_"+str(bessel_order)+"_"+ellRange+"_"+thetaRange+".table"
	file=open(Wfilename)
	W_ap=np.loadtxt(file,comments='#')
# 
	Wfilename=FolderName+"W_ap0.00_"+str(bessel_order)+"_"+ellRange+"_"+thetaRange+".table"
	file=open(Wfilename)
	W_ap0=np.loadtxt(file,comments='#')
# 
	plt.xscale("log")
	plt.plot(np.exp(W[:,0]),W[:,1],'-b')
	plt.plot(np.exp(W_ap[:,0]),W_ap[:,1],'-r')
	plt.plot(np.exp(W_ap0[:,0]),W_ap0[:,1],'-g')
	plt.axvline(x=l_min_vec[b],ls='--',color='k')
	plt.axvline(x=l_max_vec[b],ls='--',color='k')
	plt.show()



nBins=5
FolderName_outputs="example_files/outputs/"
filename=FolderName_outputs+"galaxy_cl/ell.txt"	
file=open(filename)
ell=np.loadtxt(file,comments='#')

# filename="example_files/outputs/bandpower_galaxy/l_min_vec.txt"	
# file=open(filename)
# l_min=np.loadtxt(file,comments='#')


for z1 in range(1,nBins+1):
	for z2 in range(z1,nBins+1):
		filename=FolderName_outputs+"/bandpower_galaxy/bin_"+str(z2)+"_"+str(z1)+".txt"	
		file=open(filename)
		BP_clustering=np.loadtxt(file,comments='#')
# 
		filename=FolderName_outputs+"/galaxy_cl/bin_"+str(z2)+"_"+str(z1)+".txt"	
		file=open(filename)
		Cl_clustering=np.loadtxt(file,comments='#')
		# 
		plt.xlim(100,1500)
		plt.xscale("log")
		plt.plot(ell,ell*ell*Cl_clustering,'-k')
		plt.plot(l_vec,BP_clustering,'x')
		plt.show()




for z1 in range(1,nBins+1):
	for z2 in range(z1,nBins+1):
		filename=FolderName_outputs+"/bandpower_ggl/bin_"+str(z2)+"_"+str(z1)+".txt"	
		file=open(filename)
		BP=np.loadtxt(file,comments='#')
# 
		filename=FolderName_outputs+"/galaxy_shear_cl/bin_"+str(z2)+"_"+str(z1)+".txt"	
		file=open(filename)
		Cl=np.loadtxt(file,comments='#')
		# 
		plt.xlim(100,1500)
		plt.xscale("log")
		plt.plot(ell,ell*ell*Cl,'-k')
		plt.plot(l_vec,BP,'x')
		plt.show()


for z1 in range(1,nBins+1):
	for z2 in range(z1,nBins+1):
		filename=FolderName_outputs+"/bandpower_shear_e/bin_"+str(z2)+"_"+str(z1)+".txt"	
		file=open(filename)
		BP=np.loadtxt(file,comments='#')
# 
		filename=FolderName_outputs+"/shear_cl/bin_"+str(z2)+"_"+str(z1)+".txt"	
		file=open(filename)
		Cl=np.loadtxt(file,comments='#')
		# 
		plt.xlim(100,1500)
		plt.xscale("log")
		plt.plot(ell,ell*ell*Cl,'-k')
		plt.plot(l_vec,BP,'x')
		plt.show()



# filename="table_G4_75.65.ascii"	
# file=open(filename)
# Gmu=np.loadtxt(file,comments='#')


# filename="integ_limits_noAp4_75.65.ascii"	
# file=open(filename)
# integ_lim=np.loadtxt(file,comments='#')


# # filename="integrant.ascii"	
# # file=open(filename)
# # integrant=np.loadtxt(file,comments='#')

# for b in range(nBands):
# 	filename="../../cosebis_cosmosis/integrant_Ap0.00_bin"+str(b)+"_ell_0.10_0.ascii"	
# 	file=open(filename)
# 	integrant_noap=np.loadtxt(file,comments='#')
# # 
# 	filename="../../cosebis_cosmosis/integrant_Ap0.05_bin"+str(b)+"_ell_0.10_0.ascii"	
# 	file=open(filename)
# 	integrant_ap=np.loadtxt(file,comments='#')
# # 
# 	plt.clf()
# 	plt.xscale("log")
# 	# plt.xlim(tmin*arcmin)
# 	# plt.ylim(tmax*arcmin)
# 	plt.plot(integrant_noap[:,0],integrant_noap[:,1],'-k',label="noap: bin"+str(b))
# 	plt.plot(integrant_ap[:,0],integrant_ap[:,1],'--r',label="ap: bin"+str(b))
# 	plt.show()

# # plt.clf()
# # plt.xscale("log")
# # plt.plot(integrant[:,0],integrant[:,1],'-k')
# # # plt.scatter(integ_lim[:,0],integ_lim[:,1],marker='x')
# # plt.show()

# plt.clf()
# plt.xscale("log")
# plt.plot(Gmu[:,0],Gmu[:,1],'-k')
# plt.scatter(integ_lim[:,0],integ_lim[:,1],marker='x')
# plt.show()

# plt.clf()
# bessel_order=4
# for b in range(nBands):
# 	iband=b+1
# 	gfilename="g_"+str(bessel_order)+"_"+str(iband)+"_"+thetaRange+".table"	
# 	file=open(gfilename)
# 	g=np.loadtxt(file,comments='#')
# # 
# 	gfilename="g_num_"+str(bessel_order)+"_"+str(iband)+"_"+thetaRange+".table"	
# 	file=open(gfilename)
# 	g_num=np.loadtxt(file,comments='#')
# # 
# 	plt.plot(g[:,0]/arcmin,g[:,1],'-k')
# 	plt.plot(g_num[:,0]/arcmin,g[:,1],'--r')
# 	plt.show()


# for b in range(nBands):
# 	iband=b+1
# 	Sfilename="S_"+str(iband)+"_"+str(nBands)+"_"+ellRange+".ascii"	
# 	file=open(Sfilename)
# 	S=np.loadtxt(file,comments='#')
# 	plt.clf()
# 	plt.xlim(lmin-10,lmax+100)
# 	plt.plot(S[:,0],S[:,1])
# 	plt.axvline(x=l_min_vec[b],ls='--',color='k')
# 	plt.axvline(x=l_max_vec[b],ls='--',color='k')
# 	plt.show()


