import numpy as np
import matplotlib.pylab as pl
from   matplotlib.font_manager import FontProperties
from   matplotlib.ticker import ScalarFormatter
import math, os
from matplotlib.patches import Rectangle
# from random import randint
# import random
from scipy.interpolate import interp1d
# from   array import array
from scipy import pi,sqrt,exp
# from scipy.special.orthogonal import p_roots
# from numpy.polynomial.legendre import legcompanion, legval, legder
# import numpy.linalg as la
from measure_cosebis import tplus,tminus
from argparse import ArgumentParser
from rebin import rebin




# parser = ArgumentParser(description='Take input 2pcfs files and calculate COSEBIs')
# parser.add_argument("-i", "--inputfile", dest="inputfile",
#                     help="Input file name", metavar="inputFile",required=True)

# parser.add_argument("-o", "--outputfile", dest="outputfile"
#                     ,help="output file name suffix. The outputs are En_outputfile.ascii and Bn_outputfile.ascii"
#                     , metavar="outputFile",required=True)

# parser.add_argument('-n','--nCOSEBIs', dest="nModes", type=int,default=10, nargs='?',
#          help='number of COSEBIs modes to produce')

# parser.add_argument('-s','--thetamin', dest="thetamin", type=float,default=0.5, 
#     nargs='?', help='value of thetamin')

# parser.add_argument('-l','--thetamax', dest="thetamax", type=float,default=300.0, 
#     nargs='?', help='value of thetamax')

# parser.add_argument('-c','--norm', dest="normfile", help='normalisation file name for T_plus/minus', metavar="norm",required=True)
# parser.add_argument('-r','--root', dest="rootfile", help='roots file name for T_plus/minus', metavar="root",required=True)

# # parser.add_argument('-t','--type', dest="type",default='ee',
# #  help='type of inputfile: ee: cosmic shear, ne: galaxy-galaxy lensing, nn: galaxy clustering')

# # parser.add_argument('--type_GGL', dest="type_GGL",default='tangential',
# #  help='type of GGL output: tangential or cross. Default is tangential')


# args = parser.parse_args()

# inputfile=args.inputfile
# outputfile=args.outputfile
# nModes=args.nModes
# thetamin=args.thetamin
# thetamax=args.thetamax
# normfile=arg.normfile
# rootfile=arg.rootfile

# print('input file is '+inputfile+', making COSEBIs for '+str(nModes)+' modes and theta in ['+'%.2f' %thetamin+"'" 
#     +'%.2f' %thetamax+"'], outputfiles are: En_"+outputfile+'.ascii and Bn_'+outputfile+'.ascii')



def integ_xi(xi_func,theta_edges, ntheta):
    ix=np.linspace(0,ntheta-1,ntheta)
    xip_integrated=np.zeros(len(theta_edges)-1)
    for tbin in range(len(theta_edges)-1):
        theta_in_range=np.exp(np.log(theta_edges[tbin])+(np.log(theta_edges[tbin+1])-np.log(theta_edges[tbin]))/(ntheta)*(ix+0.5))
        xip_integrated[tbin]=sum(xi_func(theta_in_range)*theta_in_range)/sum(theta_in_range)
    return xip_integrated

foldername="/Users/marika_asgary/Documents/cosmosis/MyFolders/cosebis_kids_theory_1bin/shear_xi/"
file= open(foldername+"theta.txt")
theta_rad=np.loadtxt(file,comments='#')
theta_arc=theta_rad*180.*60./np.pi

file= open(foldername+'xiplus_1_1.txt')
xip_in=np.loadtxt(file,comments='#')
xip_func=interp1d(theta_arc, xip_in, kind='cubic')

file= open(foldername+'ximinus_1_1.txt')
xim_in=np.loadtxt(file,comments='#')
xim_func=interp1d(theta_arc, xim_in, kind='cubic')

arcmin=180*60/np.pi

tmin=0.5
tmax=300.0
thetamin='%.2f' % tmin
thetamax='%.2f' % tmax
thetaRange=thetamin+'-'+thetamax  
thetamin='%.1f' % tmin
thetamax='%.0f' % tmax
thetaRangeFolder=thetamin+'-'+thetamax 
print(thetaRangeFolder)

folderName_rootnorm='/Users/marika_asgary/Documents/CosmicShear/repos/kcap/cosebis/TLogsRootsAndNorms/' 
file = open(folderName_rootnorm+'/Normalization'+'_'+thetaRange+'.table')
norm_all=np.loadtxt(file,comments='#')
# file = open()
filename=folderName_rootnorm+'/Root'+'_'+thetaRange+'.table'
roots_all = [line.strip().split() for line in open(filename)]


n=1
root_in=np.double(np.asarray(roots_all[n-1]))
root=root_in[1::]
norm_in=norm_all[n-1]
norm=norm_in[1]
tp1=tplus(tmin,tmax,n,norm,root,ntheta=10000)
tm1=tminus((tmin),(tmax),n,norm,root,tp1,ntheta=10000,nG=20)
# abs(1-Tm1_func(tm1[:,0])/tm1[:,1]).max()


n=2
root_in=np.double(np.asarray(roots_all[n-1]))
root=root_in[1::]
norm_in=norm_all[n-1]
norm=norm_in[1]
tp2=tplus(tmin,tmax,n,norm,root,ntheta=10000)
tm2=tminus(tmin,tmax,n,norm,root,tp2,ntheta=10000,nG=20*n)
# abs(1-Tm2_func(tm2[:,0])/tm2[:,1]).max()


n=5
root_in=np.double(np.asarray(roots_all[n-1]))
root=root_in[1::]
norm_in=norm_all[n-1]
norm=norm_in[1]
tp5=tplus(tmin,tmax,n,norm,root,ntheta=10000)
tm5=tminus(tmin,tmax,n,norm,root,tp5,ntheta=10000,nG=20*n)
# abs(1-Tm2_func(tm2[:,0])/tm2[:,1]).max()



n=20
root_in=np.double(np.asarray(roots_all[n-1]))
root=root_in[1::]
norm_in=norm_all[n-1]
norm=norm_in[1]
tp20=tplus(tmin,tmax,n,norm,root,ntheta=10000)
tm20=tminus(tmin,tmax,n,norm,root,tp20,ntheta=10000,nG=20*20)




folderName='/Users/marika_asgary/Documents/CosmicShear/COSEBIs/TLogs/' 

n=1
file = open(folderName+'/TpRadian'+str(n)+'_'+thetaRange+'.table')
Tp1=np.loadtxt(file,comments='#')
file = open(folderName+'/TmRadian'+str(n)+'_'+thetaRange+'.table')
Tm1=np.loadtxt(file,comments='#')
theta=np.exp(Tp1[:,0])*arcmin

Tp1_func=interp1d(theta, Tp1[:,1])
Tm1_func=interp1d(theta, Tm1[:,1])


n=2
file = open(folderName+'/TpRadian'+str(n)+'_'+thetaRange+'.table')
Tp2=np.loadtxt(file,comments='#')
file = open(folderName+'/TmRadian'+str(n)+'_'+thetaRange+'.table')
Tm2=np.loadtxt(file,comments='#')

Tp2_func=interp1d(theta, Tp2[:,1])
Tm2_func=interp1d(theta, Tm2[:,1])


n=3
file = open(folderName+'/TpRadian'+str(n)+'_'+thetaRange+'.table')
Tp3=np.loadtxt(file,comments='#')
file = open(folderName+'/TmRadian'+str(n)+'_'+thetaRange+'.table')
Tm3=np.loadtxt(file,comments='#')

Tp3_func=interp1d(theta, Tp3[:,1])
Tm3_func=interp1d(theta, Tm3[:,1])


n=4
file = open(folderName+'/TpRadian'+str(n)+'_'+thetaRange+'.table')
Tp4=np.loadtxt(file,comments='#')
file = open(folderName+'/TmRadian'+str(n)+'_'+thetaRange+'.table')
Tm4=np.loadtxt(file,comments='#')



Tp4_func=interp1d(theta, Tp4[:,1])
Tm4_func=interp1d(theta, Tm4[:,1])


n=5
file = open(folderName+'/TpRadian'+str(n)+'_'+thetaRange+'.table')
Tp5=np.loadtxt(file,comments='#')
file = open(folderName+'/TmRadian'+str(n)+'_'+thetaRange+'.table')
Tm5=np.loadtxt(file,comments='#')

Tp5_func=interp1d(theta, Tp5[:,1])
Tm5_func=interp1d(theta, Tm5[:,1])


n=10
file = open(folderName+'/TpRadian'+str(n)+'_'+thetaRange+'.table')
Tp10=np.loadtxt(file,comments='#')
file = open(folderName+'/TmRadian'+str(n)+'_'+thetaRange+'.table')
Tm10=np.loadtxt(file,comments='#')

Tp10_func=interp1d(theta, Tp10[:,1])
Tm10_func=interp1d(theta, Tm10[:,1])


n=20
file = open(folderName+'/TpRadian'+str(n)+'_'+thetaRange+'.table')
Tp20=np.loadtxt(file,comments='#')
file = open(folderName+'/TmRadian'+str(n)+'_'+thetaRange+'.table')
Tm20=np.loadtxt(file,comments='#')

Tp20_func=interp1d(theta, Tp20[:,1])
Tm20_func=interp1d(theta, Tm20[:,1])


n=7
file = open(folderName+'/TpRadian'+str(n)+'_'+thetaRange+'.table')
Tp7=np.loadtxt(file,comments='#')
file = open(folderName+'/TmRadian'+str(n)+'_'+thetaRange+'.table')
Tm7=np.loadtxt(file,comments='#')

Tp7_func=interp1d(theta, Tp7[:,1])
Tm7_func=interp1d(theta, Tm7[:,1])


iterations=20
Eplus_vec=np.zeros(iterations)
Eminus_vec=np.zeros(iterations)
Eplus_vec_py=np.zeros(iterations)
Eminus_vec_py=np.zeros(iterations)

Eplus_vec_lin=np.zeros(iterations)
Eminus_vec_lin=np.zeros(iterations)
nbins_vec=np.linspace(100,1000,iterations)

tp2_func=interp1d(tp2[:,0], tp2[:,1])
tm2_func=interp1d(tm2[:,0], tm2[:,1])

tp1_func=interp1d(tp1[:,0], tp1[:,1])
tm1_func=interp1d(tm1[:,0], tm1[:,1])

tp5_func=interp1d(tp5[:,0], tp5[:,1])
tm5_func=interp1d(tm5[:,0], tm5[:,1])


tp20_func=interp1d(tp20[:,0], tp20[:,1])
tm20_func=interp1d(tm20[:,0], tm20[:,1])


for ibins in range(iterations):
    nbins=nbins_vec[ibins]
    ix=np.linspace(0,nbins-1,nbins)
    theta_mid=np.exp(np.log(tmin)+(np.log(tmax)-np.log(tmin))/(nbins)*(ix+0.5))
    theta_edges=np.logspace(np.log10(tmin),np.log10(tmax),nbins+1)
# 
    theta_low=theta_edges[0:-1]
    theta_high=theta_edges[1::]
    delta_theta=theta_high-theta_low
#
    integ_plus=Tp5_func(theta_mid)*theta_mid*integ_xi(xip_func,theta_edges,100)
    integ_minus=Tm5_func(theta_mid)*theta_mid*integ_xi(xim_func,theta_edges,100)
# 
    integ_plus_py=tp5_func(theta_mid)*theta_mid*integ_xi(xip_func,theta_edges,100)
    integ_minus_py=tm5_func(theta_mid)*theta_mid*integ_xi(xim_func,theta_edges,100)
# 
# 
    Eplus_vec[ibins]=sum(integ_plus*delta_theta)
    Eminus_vec[ibins]=sum(integ_minus*delta_theta)
# 
    Eplus_vec_py[ibins]=sum(integ_plus_py*delta_theta)
    Eminus_vec_py[ibins]=sum(integ_minus_py*delta_theta)
# 
    theta_edges_lin=np.linspace(tmin,tmax,nbins+1)
    theta_lin_low=theta_edges_lin[0:-1]
    theta_lin_high=theta_edges_lin[1::]
    theta_lin=(theta_lin_high+theta_lin_low)/2.
    delta_theta_lin=(theta_lin_high-theta_lin_low)
# 
    integ_plus_lin=Tp5_func(theta_lin)*theta_lin*integ_xi(xip_func,theta_edges_lin,100)
    integ_minus_lin=Tm5_func(theta_lin)*theta_lin*integ_xi(xim_func,theta_edges_lin,100)
# 
# 
    Eplus_vec_lin[ibins]=sum(integ_plus_lin*delta_theta_lin)
    Eminus_vec_lin[ibins]=sum(integ_minus_lin*delta_theta_lin)


pl.clf()
pl.plot(nbins_vec,np.zeros(len(nbins_vec)),'-k')
# pl.plot(nbins_vec,Eplus_vec_py/Eplus_vec-1,'--r')
# pl.plot(nbins_vec,Eminus_vec_py/Eminus_vec-1,'--g')
# pl.plot(nbins_vec,0.5*(Eplus_vec+Eminus_vec)/(0.5*(Eplus_vec[-1]+Eminus_vec[-1]))-1,'--r')
pl.plot(nbins_vec,0.5*(Eplus_vec_py+Eminus_vec_py)/(0.5*(Eplus_vec[-1]+Eminus_vec[-1]))-1,'--g',label='En convergence')
pl.plot(nbins_vec,Eplus_vec-Eminus_vec,'-r',label='c++')
pl.plot(nbins_vec,Eplus_vec_py-Eminus_vec_py,'--b',label='python')
# pl.plot(nbins_vec,0.5*(Eplus_vec_lin+Eminus_vec_lin)/(0.5*(Eplus_vec[-1]+Eminus_vec[-1]))-1,'b')
# pl.plot(nbins_vec,Eminus_vec,'b')
pl.legend(loc='best')
pl.show()

# pl.clf()
# # pl.plot(nbins_vec,0.5*(Eplus_vec+Eminus_vec)/(0.5*(Eplus_vec[-1]+Eminus_vec[-1]))-1,'r')
# pl.plot(nbins_vec,0.5*(Eplus_vec+Eminus_vec)/(0.5*(Eplus_vec[-1]+Eminus_vec[-1]))-1,'r')
# pl.plot(nbins_vec,0.5*(Eplus_vec_lin+Eminus_vec_lin)/(0.5*(Eplus_vec[-1]+Eminus_vec[-1]))-1,'b')
# # pl.plot(nbins_vec,Eminus_vec,'b')
# pl.show()


pl.rcParams.update(params)
pl.clf()
pl.figure(1)
pl.subplots_adjust(wspace=0,hspace=0)
yprops=dict(rotation=90, horizontalalignment='center',verticalalignment='center',x=10,fontsize=13,labelpad=10)
leg1=Rectangle((0,0),0,0,alpha=0.0)

ax=pl.subplot(1,1,1)
ax.set_xlim(0.49,100)
ax.set_xscale('log')
pl.xlabel(r'$\theta(arcmin)$')
pl.ylabel(r'$T_{+n}(\theta)$',**yprops)
ax.plot(theta,Tp5[:,1],':r',lw=2.0, label=r'$n=5$')
ax.plot(theta_mid,Tp5_func(theta_mid),lw=2.0, label=r'interpolated')

lg=pl.legend(loc='lower right')
pl.gca().add_artist(lg)
lg1 = pl.legend([leg1],[r'['+thetamin+r"',"+thetamax+r"']"],
            bbox_to_anchor=(0.2, 0.2),handlelength=0,
            borderpad=0.5,labelspacing=0.4,numpoints=1,prop={'size':13})
lg1.draw_frame(False)
pl.show()



pl.rcParams.update(params)
pl.clf()
pl.figure(1)
pl.subplots_adjust(wspace=0,hspace=0)
yprops=dict(rotation=90, horizontalalignment='center',verticalalignment='center',x=10,fontsize=13,labelpad=10)
leg1=Rectangle((0,0),0,0,alpha=0.0)


nbins=100
ix=np.linspace(0,nbins-1,nbins)
theta_mid=np.exp(np.log(tmin)+(np.log(tmax)-np.log(tmin))/(nbins)*(ix+0.5))
theta_edges=np.logspace(np.log10(tmin),np.log10(tmax),nbins+1)

integrand1=Tp1_func(theta_mid)*theta_mid*integ_xi(xip_func,theta_edges,100)
integrand2=Tp2_func(theta_mid)*theta_mid*integ_xi(xip_func,theta_edges,100)
integrand3=Tp3_func(theta_mid)*theta_mid*integ_xi(xip_func,theta_edges,100)
integrand4=Tp4_func(theta_mid)*theta_mid*integ_xi(xip_func,theta_edges,100)
integrand5=Tp5_func(theta_mid)*theta_mid*integ_xi(xip_func,theta_edges,100)

integ_lin5=Tp5_func(theta_lin)*theta_lin*integ_xi(xip_func,theta_edges_lin,100)

ax=pl.subplot(1,1,1)
ax.set_xlim(0.49,300)
ax.set_xscale('log')
pl.xlabel(r'$\theta(arcmin)$')
pl.ylabel(r'$T_{+n}(\theta)\xi_+(\theta)\theta$',**yprops)
ax.plot(theta,Tp1[:,1]*theta*xip_func(theta),'-r',lw=2.0, label=r'$n=1$')
ax.plot(theta_mid,integrand1,'m',marker='o',markerfacecolor='none',lw=0.5,mew=1, label=r'integ1')
ax.plot(theta,Tp2[:,1]*theta*xip_func(theta),'-b',lw=2.0, label=r'$n=2$')
ax.plot(theta_mid,integrand2,'c',marker='o',markerfacecolor='none',lw=0.5,mew=1, label=r'integ2')
ax.plot(theta,Tp3[:,1]*theta*xip_func(theta),'-y',lw=2.0, label=r'$n=3$')
ax.plot(theta_mid,integrand3,'g',marker='o',markerfacecolor='none',lw=0.5,mew=1, label=r'integ3')
ax.plot(theta,Tp4[:,1]*theta*xip_func(theta),'-k',lw=2.0, label=r'$n=4$')
ax.plot(theta_mid,integrand4,'r',marker='o',markerfacecolor='none',lw=0.5,mew=1, label=r'integ4')
ax.plot(theta,Tp5[:,1]*theta*xip_func(theta),'-g',lw=2.0, label=r'$n=5$')
ax.plot(theta_mid,integrand5,'b',marker='o',markerfacecolor='none',lw=0.5,mew=1, label=r'integ5')


lg=pl.legend(loc='lower right')
pl.gca().add_artist(lg)
lg1 = pl.legend([leg1],[r'['+thetamin+r"',"+thetamax+r"']"],
            bbox_to_anchor=(0.2, 0.2),handlelength=0,
            borderpad=0.5,labelspacing=0.4,numpoints=1,prop={'size':13})
lg1.draw_frame(False)
pl.show()



n=10
root_in=np.double(np.asarray(roots_all[n-1]))
root=root_in[1::]
norm_in=norm_all[n-1]
norm=norm_in[1]
tp5=tplus(tmin,tmax,n,norm,root,ntheta=10000)
tm5=tminus(tmin,tmax,n,norm,root,tp5,ntheta=10000,nG=20*n)
tm5_quad=tminus_quad(tmin,tmax,n,norm,root,tp5,ntheta=10000)


tp5_func=interp1d(tp5[:,0], tp5[:,1])
tm5_func=interp1d(tm5[:,0], tm5[:,1])

pl.clf()
pl.xlim(tmin,tmax)
pl.xscale('log')
pl.xlabel(r'$\theta(arcmin)$')
pl.ylabel(r'$T_{-n}(\theta)$')
# pl.plot(tm5[:,0],tm5_func(tm5[:,0]),lw=2.0, label=r'python')
pl.plot(tm5_quad[:,0],tm5_quad[:,1],'--k',lw=2.0, label=r'quad')
pl.plot(tm5[:,0],Tm10_func(tm5[:,0]),':r',lw=2.0, label=r'$n=$'+str(5))
pl.legend(loc='best')
pl.show()


# pl.savefig('plots/Tpn'+thetaRangeFolder+'.pdf',bbox_inches='tight')
