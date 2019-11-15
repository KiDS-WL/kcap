import numpy as np
# import matplotlib.pylab as pl
# from   matplotlib.font_manager import FontProperties
# from   matplotlib.ticker import ScalarFormatter
import math, os
# from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d
from scipy import pi,sqrt,exp
from measure_cosebis import tminus_quad, tplus,tminus
from argparse import ArgumentParser
from rebin import rebin



# example run:
# tmin=0.50
# tmax=100.00
# python run_measure_cosebis_cats2stats.py -i xi_nBins_1_Bin1_Bin1 -o nBins_1_Bin1_Bin1 --norm ./TLogsRootsAndNorms/Normalization_${tmin}-${tmax}.table -r ./TLogsRootsAndNorms/Root_${tmin}-${tmax}.table -b lin --thetamin ${tmin} --thetamax ${tmax} -n 20

parser = ArgumentParser(description='Take input 2pcfs files and calculate COSEBIs')
parser.add_argument("-i", "--inputfile", dest="inputfile",
    help="Full Input file name", metavar="inputFile",required=True)

parser.add_argument('-t','--theta_col', dest="theta_col", type=int,default=0, nargs='?',
         help='column for theta, default is 0')

parser.add_argument('-p','--xip_col', dest="xip_col", type=int,default=1, nargs='?',
         help='column for xi_plus, default is 1')

parser.add_argument('-m','--xim_col', dest="xim_col", type=int,default=2, nargs='?',
         help='column for xi_minus, default is 2')


parser.add_argument('--cfoldername', dest="cfoldername", 
    help='full name and address of the folder for En/Bn files, default is cosebis_results',default="./cosebis_results",required=False)

parser.add_argument("-o", "--outputfile", dest="outputfile"
                    ,help="output file name suffix. The outputs are cfoldername/En_${outputfile}.ascii and cfoldername/Bn_${outputfile}.ascii"
                    ,metavar="outputFile",required=True)


parser.add_argument('-b','--binning', dest="binning", help='log or lin binning, default is log',default="log",required=False)

parser.add_argument('-n','--nCOSEBIs', dest="nModes", type=int,default=10, nargs='?',
         help='number of COSEBIs modes to produce, default is 10')

parser.add_argument('-s','--thetamin', dest="thetamin", type=float,default=0.5, 
    nargs='?', help='value of thetamin in arcmins')

parser.add_argument('-l','--thetamax', dest="thetamax", type=float,default=300.0, 
    nargs='?', help='value of thetamax, in arcmins')

parser.add_argument('--tfoldername', dest="tfoldername", help='name and full address of the folder for Tplus Tminus files, will make it if it does not exist',default="Tplus_minus",required=False)
parser.add_argument('--tplusfile', dest="tplusfile", help='name of Tplus file, will look for it before running the code',default="Tplus",required=False)
parser.add_argument('--tminusfile', dest="tminusfile", help='name of Tplus file, will look for it before running the code',default="Tminus",required=False)

parser.add_argument('-c','--norm', dest="normfile", help='normalisation file name and address for T_plus/minus', metavar="norm",required=True)
parser.add_argument('-r','--root', dest="rootfile", help='roots file name and address for T_plus/minus', metavar="root",required=True)


args = parser.parse_args()

inputfile=args.inputfile
theta_col=args.theta_col
xip_col=args.xip_col
xim_col=args.xim_col
outputfile=args.outputfile
nModes=args.nModes
thetamin=args.thetamin
thetamax=args.thetamax
normfile=args.normfile
rootfile=args.rootfile
tplusfile=args.tplusfile
tminusfile=args.tminusfile
tfoldername=args.tfoldername
cfoldername=args.cfoldername
binning=args.binning

print('input file is '+inputfile+', making COSEBIs for '+str(nModes)+' modes and theta in ['+'%.2f' %thetamin+"'," 
    +'%.2f' %thetamax+"'], outputfiles are: "+cfoldername+"/En_"+outputfile+'.ascii and '+cfoldername+'/Bn_'+outputfile+'.ascii')


file=open(inputfile)
xipm_in=np.loadtxt(file,comments='#')
theta=xipm_in[:,theta_col]
xip=xipm_in[:,xip_col]
xim=xipm_in[:,xim_col]
np_gals = xipm_in[:,7]


if(binning=='log'):
    good_args=np.squeeze(np.argwhere((theta>thetamin) & (theta<thetamax)))
    nbins_within_range=len(theta[good_args])
    ix=np.linspace(0,nbins_within_range-1,nbins_within_range)
    theta_mid=np.exp(np.log(thetamin)+(np.log(thetamax)-np.log(thetamin))/(nbins_within_range)*(ix+0.5))
# 
    theta_edges=np.logspace(np.log10(thetamin),np.log10(thetamax),nbins_within_range+1)
    # 
    theta_low=theta_edges[0:-1]
    theta_high=theta_edges[1::]
    delta_theta=theta_high-theta_low
# 
    #check if the mid points are close enough
    if((abs(theta_mid/theta[good_args]-1)<(delta_theta/10.)).all()):
        print("I'm happy")
    else:
        print("Not happy with input thetas, exiting now ...")
        exit()
        # 
elif(binning=='lin'):
    good_args=np.squeeze(np.argwhere((theta>thetamin) & (theta<thetamax)))
    nbins_within_range=len(theta[good_args])
    delta_theta=np.zeros(nbins_within_range)
    delta_theta[1::]=(theta[1]-theta[0])
    delta_theta[0]=(theta[1]-theta[0])/2.
    delta_theta[-1]=(theta[1]-theta[0])/2.
    theta_mid=np.linspace(thetamin+delta_theta[0],thetamax-delta_theta[-1],nbins_within_range)
    # 
elif(binning=='lin_to_log'):
    nbins=300
    xip_binned=rebin(theta,xip,np_gals,thetamin,thetamax,nbins)
    xim_binned=rebin(theta,xim,np_gals,thetamin,thetamax,nbins)
    ix=np.linspace(0,nbins-1,nbins)
    theta_mid=np.exp(np.log(thetamin)+(np.log(thetamax)-np.log(thetamin))/(nbins)*(ix+0.5))
    xip=xip_binned[:,2]
    xim=xim_binned[:,2]
    theta=xip_binned[:,0]
    good_args=np.squeeze(np.argwhere((theta>thetamin) & (theta<thetamax)))
    theta_edges=np.logspace(np.log10(thetamin),np.log10(thetamax),nbins+1)
    # 
    theta_low=theta_edges[0:-1]
    theta_high=theta_edges[1::]
    delta_theta=theta_high-theta_low
    # 
else:
    print('not a recognised binning scheme, exiting now')
    exit()

arcmin=180*60/np.pi

tmin='%.2f' % thetamin
tmax='%.2f' % thetamax
thetaRange=tmin+'-'+tmax  


if not os.path.exists(tfoldername):
    os.makedirs(tfoldername)

if not os.path.exists(cfoldername):
    os.makedirs(cfoldername)

En=np.zeros(nModes)
Bn=np.zeros(nModes)

for n in range(1,nModes+1):
    # 
    TplusFileName=tfoldername+"/"+tplusfile+str(n)+"_"+thetaRange+".table"
    TminusFileName=tfoldername+"/"+tminusfile+str(n)+"_"+thetaRange+".table"
    if(os.path.isfile(TplusFileName)):
        file = open(TplusFileName)
        tp=np.loadtxt(file,comments='#')
    else:
        file = open(normfile)
        norm_all=np.loadtxt(file,comments='#')
        norm_in=norm_all[n-1]
        norm=norm_in[1]
        # 
        roots_all = [line.strip().split() for line in open(rootfile)]
        root_in=np.double(np.asarray(roots_all[n-1]))
        root=root_in[1::]
        # 
        tp=tplus(thetamin,thetamax,n,norm,root,ntheta=10000)
        np.savetxt(TplusFileName,tp)
# 
    if(os.path.isfile(TminusFileName)):
        file = open(TminusFileName)
        tm=np.loadtxt(file,comments='#')
    else:
        file = open(normfile)
        norm_all=np.loadtxt(file,comments='#')
        norm_in=norm_all[n-1]
        norm=norm_in[1]
    # 
        roots_all = [line.strip().split() for line in open(rootfile)]
        root_in=np.double(np.asarray(roots_all[n-1]))
        root=root_in[1::]
        # 
        tm=tminus(thetamin,thetamax,n,norm,root,tp,ntheta=10000)
        np.savetxt(TminusFileName,tm)

    tp_func=interp1d(tp[:,0], tp[:,1])
    tm_func=interp1d(tm[:,0], tm[:,1])

    ##plots the T_p/m 
    # pl.clf()
    # # pl.xlim(thetamin,thetamax)
    # pl.xscale('log')
    # pl.xlabel(r'$\theta(arcmin)$')
    # pl.ylabel(r'$T_{+n}(\theta)$')
    # pl.plot(tp[:,0],tp_func(tp[:,0]),lw=2.0, label=r'python')
    # pl.plot(theta_t,Tp[:,1],':r',lw=2.0, label=r'$n=$'+str(n))
    # # pl.plot(tp[:,0],tp[:,1],lw=2.0, label=r'python')
    # pl.legend(loc='best')
    # pl.show()
    # # 
    # pl.clf()
    # # pl.xlim(thetamin,thetamax)
    # pl.xscale('log')
    # pl.xlabel(r'$\theta(arcmin)$')
    # pl.ylabel(r'$T_{-n}(\theta)$')
    # pl.plot(tm[:,0],tm_func(tm[:,0]),lw=2.0, label=r'python')
    # pl.plot(tm[:,0],Tm_func(tm[:,0]),':r',lw=2.0, label=r'$n=$'+str(n))
    # # pl.plot(tm[:,0],tm[:,1],lw=2.0, label=r'python')
    # pl.legend(loc='best')
    # pl.show()
    # 
    integ_plus=tp_func(theta_mid)*theta_mid*xip[good_args]
    integ_minus=tm_func(theta_mid)*theta_mid*xim[good_args]
    # 
    Integral_plus=sum(integ_plus*delta_theta)
    Integral_minus=sum(integ_minus*delta_theta)
    En[n-1]=0.5*(Integral_plus+Integral_minus)/arcmin/arcmin
    Bn[n-1]=0.5*(Integral_plus-Integral_minus)/arcmin/arcmin

EnfileName=cfoldername+"/En_"+outputfile+".ascii"
BnfileName=cfoldername+"/Bn_"+outputfile+".ascii"
np.savetxt(EnfileName,En)
np.savetxt(BnfileName,Bn)



