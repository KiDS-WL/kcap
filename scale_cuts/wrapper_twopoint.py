

    ##############################
    ##  wrapper_twopoint.py	##
    ##  Chieh-An Lin		##
    ##  Version 2019.11.11	##
    ##############################


import numpy as np
import scipy as sp
import astropy.io.fits as fits
import twopoint


###############################################################################
## Functions related to label convention

class LabelConvention:
  
    def __init__(self, w='wTh', gamma_t='gT', gamma_x='gX', xi_p='xiP', xi_m='xiM', 
                  P_nn='Pnn', P_ne_E='PneE', P_ne_B='PneB', P_ee_E='PeeE', P_ee_B='PeeB', 
                  E_n='En', B_n='Bn'):
    
        self.prefix  = 'twoPoint'
        self.lens    = 'NZ_LENS'
        self.source  = 'NZ_SOURCE'
        
        self.w       = 'wTh'
        self.gamma_t = 'gT'
        self.gamma_x = 'gX'
        self.xi_p    = 'xiP'
        self.xi_m    = 'xiM'
        self.P_nn    = 'Pnn'
        self.P_ne_E  = 'PneE'
        self.P_ne_B  = 'PneB'
        self.P_ee_E  = 'PeeE'
        self.P_ee_B  = 'PeeB'
        self.E_n     = 'En'
        self.B_n     = 'Bn'
        
        self.labDict = {
            self.w:        w,
            self.gamma_t:  gamma_t,
            self.gamma_x:  gamma_x,
            self.xi_p:     xi_p,
            self.xi_m:     xi_m,
            self.P_nn:     P_nn,
            self.P_ne_E:   P_ne_E,
            self.P_ne_B:   P_ne_B,
            self.P_ee_E:   P_ee_E,
            self.P_ee_B:   P_ee_B,
            self.E_n:      E_n,
            self.B_n:      B_n,
            'all':         'all'
        }
      
        self.revLabDict = {}
        for key, value in self.labDict.items():
            self.revLabDict[value] = key
        
        tpType1 = twopoint.Types.galaxy_position_real       ## GPR
        tpType2 = twopoint.Types.galaxy_shear_plus_real     ## G+R
        tpType3 = twopoint.Types.galaxy_shear_minus_real    ## G-R
        tpType4 = twopoint.Types.galaxy_position_fourier    ## GPF
        tpType5 = twopoint.Types.galaxy_shear_emode_fourier ## GEF
        tpType6 = twopoint.Types.galaxy_shear_bmode_fourier ## GBF
        
        self.kernelTypeDict = { ## Don't touch the order of this list
            w:       [self.lens,   self.lens,   tpType1, tpType1, 'arcmin'], 
            gamma_t: [self.lens,   self.source, tpType1, tpType2, 'arcmin'], 
            gamma_x: [self.lens,   self.source, tpType1, tpType3, 'arcmin'], 
            xi_p:    [self.source, self.source, tpType2, tpType2, 'arcmin'],  
            xi_m:    [self.source, self.source, tpType3, tpType3, 'arcmin'],
            P_nn:    [self.lens,   self.lens,   tpType4, tpType4, None], 
            P_ne_E:  [self.lens,   self.source, tpType4, tpType5, None], 
            P_ne_B:  [self.lens,   self.source, tpType4, tpType6, None], 
            P_ee_E:  [self.source, self.source, tpType5, tpType5, None], 
            P_ee_B:  [self.source, self.source, tpType6, tpType6, None], 
            E_n:     [self.source, self.source, tpType5, tpType6, None], 
            B_n:     [self.source, self.source, tpType6, tpType5, None]
        }
        return
    
    def defaultToCustomStatsTag(self, statsTag):
        statsList   = statsTag.split('+')
        statsList_c = [self.labDict[stats] for stats in statsList]
        statsTag_c  = '+'.join(statsList_c)
        return statsTag_c
    
    def customToDefaultStatsTag(self, statsTag_c):
        statsList_c = statsTag_c.split('+')
        statsList   = [self.revLabDict[stats_c] for stats_c in statsList_c]
        statsTag    = '+'.join(statsList)
        return statsTag
    
    def makeScaleCutsArgs(self, scDict):
        statsList = None
        cutCross = False
        statsTag_tomoInd_tomoInd_list = []
        statsTag_binIndList_dict = {}
        statsTag_tomoInd1_tomoInd2__angMin_angMax_dict = {}
        statsTag__angMin_angMax_dict = {}
      
        for key, value in scDict.items():
            keySplit = key.split('_')
            
            if key == 'use_stats':
                statsList = value.split()
            
            elif key == 'cut_cross':
                cutCross = bool(value)
            
            elif 'cut_pair' in key:
                statsTag = self.labDict[keySplit[2]]
                for pair in value.split():
                    pair = pair.split('+')
                    tomoInd1 = int(pair[0])
                    tomoInd2 = int(pair[1])
                    statsTag_tomoInd_tomoInd_list.append((statsTag, tomoInd1, tomoInd2))
            
            elif 'cut_bin' in key:
                statsTag   = self.labDict[keySplit[2]]
                binIndList = [int(ind) for ind in value.split()]
                statsTag_binIndList_dict[statsTag] = binIndList
            
            elif 'keep_ang' in key and len(keySplit) > 4:
                statsTag = self.labDict[keySplit[2]]
                tomoInd1 = int(keySplit[3])
                tomoInd2 = int(keySplit[4])
                value    = value.split()
                angMin   = float(value[0])
                angMax   = float(value[1])
                statsTag_tomoInd1_tomoInd2__angMin_angMax_dict[(statsTag, tomoInd1, tomoInd2)] = (angMin, angMax)
            
            elif 'keep_ang' in key and len(keySplit) < 4:
                statsTag = self.labDict[keySplit[2]]
                value    = value.split()
                angMin   = float(value[0])
                angMax   = float(value[1])
                statsTag__angMin_angMax_dict[statsTag] = (angMin, angMax)
            
        scArgs = cutCross, statsTag_tomoInd_tomoInd_list, statsTag_binIndList_dict, \
            statsTag_tomoInd1_tomoInd2__angMin_angMax_dict, statsTag__angMin_angMax_dict
        return statsList, scArgs

###############################################################################
## Functions related to SpectrumMeasurement builder

class SpectrumBuilder():
  
    def __init__(self):
        self.tIL1    = []
        self.tIL2    = []
        self.aIL     = []
        self.angList = []
        self.valList = []
        return
    
    def addTomo(self, tomoInd1, tomoInd2, angle, value):
        N_ang  = len(angle)
        angInd = np.arange(N_ang, dtype=int)
        self.tIL1.append([tomoInd1+1] * N_ang)
        self.tIL2.append([tomoInd2+1] * N_ang)
        self.aIL.append(angInd + 1)
        self.angList.append(angle)
        self.valList.append(value)
        return
    
    def makeSpectrum(self, name, types, angle_unit, kernels=(None, None)):
        self.tIL1    = np.concatenate(self.tIL1)
        self.tIL2    = np.concatenate(self.tIL2)
        self.aIL     = np.concatenate(self.aIL)
        self.angList = np.concatenate(self.angList)
        self.valList = np.concatenate(self.valList)
        spec = twopoint.SpectrumMeasurement(name, (self.tIL1, self.tIL2), types, kernels, 'SAMPLE', self.aIL, self.valList, angle=self.angList, angle_unit=angle_unit)
        return spec

###############################################################################
## Functions related to TwoPointFile wrapper

class TwoPointWrapper(twopoint.TwoPointFile):
  
    def __init__(self, spectra, kernels, windows, covmat_info):
        super().__init__(spectra, kernels, windows, covmat_info)
        return
    
    @classmethod
    def from_spectra(cls, spectra, kernels=None, covmat_info=None):
        TP = cls(spectra, kernels, 'SAMPLE', covmat_info) ## windows = 'SAMPLE'
        return TP
    
    @classmethod
    def from_fits(cls, filename, covmat_name=None):
        TP = super(TwoPointWrapper, cls).from_fits(filename, covmat_name=covmat_name)
        return TP
    
    def makeMeanVector(self):
        theory = []
        for spectrum in self.spectra:
            theory.append(spectrum.value)
        theory = np.concatenate(theory)
        return theory
      
    def replaceMeanVector(self, data):
        idx = []
        for spectrum in self.spectra:
            idx.append(len(spectrum.value))
        
        if sum(idx) != len(data):
            raise ValueError(f"Size of provided data vector incompatible with spectra: {len(data)} vs {sum(idx)}.")

        idx.insert(0, 0)
        idx = np.cumsum(idx)

        for i, spectrum in enumerate(self.spectra):
            spectrum.value = data[idx[i]:idx[i+1]]
        
    def cutScales(self, cutCross=False, statsTag_tomoInd_tomoInd_list=[], statsTag_binIndList_dict={}):
        if cutCross:
            self.mask_cross() ## Cut cross pair bins, but not cross-corr between auto pair bins; e.g. C_12_12 masked, C_11_22 preserved
        
        if len(statsTag_tomoInd_tomoInd_list) > 0:
            self.mask_scales(bin_cuts=statsTag_tomoInd_tomoInd_list) ## Cut the whole tomo pair bin
        
        for statsTag, binIndList in statsTag_binIndList_dict.items():
            self.mask_indices(statsTag, binIndList) ## Cut by bin indices, 0-indexing
        return
    
    def keepScales(self, statsTag_tomoInd1_tomoInd2__angMin_angMax_dict={}, statsTag__angMin_angMax_dict={}):
        if len(statsTag_tomoInd1_tomoInd2__angMin_angMax_dict) > 0:
            self.mask_scales(cuts=statsTag_tomoInd1_tomoInd2__angMin_angMax_dict) ## Keep an angular range for specific tomo bins
          
        for statsTag, angMin_angMax in statsTag__angMin_angMax_dict.items():
            self.mask_scale(statsTag, min_scale=angMin_angMax[0], max_scale=angMin_angMax[1]) ## Keep an angular range for a specific tomo bin
        return

###############################################################################
## Functions related to mean & covariance

class TwoPointBuilder:
  
    def __init__(self,
                  nbTomoN=2, 
                  nbTomoG=5, 
                  nOfZNameList=None, 
                  nGalList=None, 
                  sigmaEpsList=None, 
                  N_theta=9, 
                  theta_min=0.5, 
                  theta_max=300,
                  N_ell=8, 
                  ell_min=100, 
                  ell_max=1500,
                  nbModes=5,
                  prefix_Flinc=None,
                  prefix_CosmoSIS=None,
                  verbose=True):
    
        ## n(z), n_gal, & sigma_eff input - otis
        self.nOfZNameList = [
            'data/mockFootprint/milvus/MFP_selection/nOfZ_hist_BOSSA_tomo0.dat',
            'data/mockFootprint/milvus/MFP_selection/nOfZ_hist_BOSSA_tomo1.dat',
            'data/mockFootprint/milvus/MFP_selection/nOfZ_hist_KiDSVDA_tomo0.dat',
            'data/mockFootprint/milvus/MFP_selection/nOfZ_hist_KiDSVDA_tomo1.dat',
            'data/mockFootprint/milvus/MFP_selection/nOfZ_hist_KiDSVDA_tomo2.dat',
            'data/mockFootprint/milvus/MFP_selection/nOfZ_hist_KiDSVDA_tomo3.dat',
            'data/mockFootprint/milvus/MFP_selection/nOfZ_hist_KiDSVDA_tomo4.dat'
        ] if nOfZNameList is None else nOfZNameList
        self.nGalList     = [0.016702, 0.016725,   0.85651964, 1.5625674, 2.2357926, 1.5223054, 1.384634]   if nGalList is None else nGalList
        self.sigmaEpsList = ([0.284605] * nbTomoG)   if sigmaEpsList is None else sigmaEpsList
        
        ## File prefixes
        self.prefix_Flinc    = 'data/mockFootprint/otis/MFP_combDataMat/'   if prefix_Flinc is None else prefix_Flinc
        self.prefix_CosmoSIS = 'data/mockFootprint/otis/cosmosis/'   if prefix_CosmoSIS is None else prefix_CosmoSIS
        
        ## Customize the above for your own inference; but don't touch the below
        ########################################################################
        
        nbTomo = nbTomoN + nbTomoG
        assert nbTomo == len(self.nOfZNameList),  "Bad length of nOfZNameList"
        assert nbTomo == len(self.nGalList), "Bad length of nGalList"
      
        if len(self.sigmaEpsList) == nbTomoG:
            pass
        elif len(self.sigmaEpsList) == nbTomo:
            self.sigmaEpsList = self.sigmaEpsList[nbTomoN:]
        else:
            raise ValueError('Bad length of sigmaEpsList')
        
        ## Tomographic bins
        self.nbTomoN = nbTomoN
        self.nbTomoG = nbTomoG
        
        ## Define angular bin parameters
        self.N_theta   = N_theta
        self.theta_min = theta_min
        self.theta_max = theta_max
        self.N_ell     = N_ell
        self.ell_min   = ell_min
        self.ell_max   = ell_max
        self.nbModes   = nbModes
        
        self.verbose    = verbose
        self.kernelList = None
        
        self.setAngBins()
        self.setNbPairs()
        self.setKernels()
        return
    
    ## Initialize
    def setAngBins(self):
        bAndC_theta = np.logspace(np.log10(self.theta_min), np.log10(self.theta_max), 2*self.N_theta+1)
        self.theta     = bAndC_theta[1::2] ## [arcmin]
        self.bin_theta = bAndC_theta[0::2]
        bAndC_ell = np.logspace(np.log10(self.ell_min), np.log10(self.ell_max), 2*self.N_ell+1)
        self.ell       = bAndC_ell[1::2]
        self.bin_ell   = bAndC_ell[0::2]
        self.nArr      = np.arange(self.nbModes) + 1.0
        return
    
    def setNbPairs(self):
        self.__nbPairsNN  = self.nbTomoN * (self.nbTomoN+1) // 2
        self.__nbPairsNG  = self.nbTomoN * self.nbTomoG
        self.__nbPairsGG  = self.nbTomoG * (self.nbTomoG+1) // 2
        self.__pairListNN = [(i, j) for i in range(self.nbTomoN) for j in range(i, self.nbTomoN)]
        self.__pairListNG = [(i, j) for i in range(self.nbTomoN) for j in range(self.nbTomoG)]
        self.__pairListGG = [(i, j) for i in range(self.nbTomoG) for j in range(i, self.nbTomoG)]
        return
    
    ## Kernel
    def _makeKernel(self, name, nOfZNameList, nGalList, sigmaEpsList):
        if len(nOfZNameList) == 0:
            return None
        
        nOfZName = nOfZNameList[0]
        if nOfZName[-4:] == '.npy':
            nOfZList = [np.load(nOfZName) for nOfZName in nOfZNameList]
        elif nOfZName[-4:] == '.fit' or nOfZName[-5:] == '.fits':
            nOfZList = [fits.getdata(nOfZName, 1) for nOfZName in nOfZNameList]
            nOfZList = [[data.field(0), data.field(1)] for data in nOfZList]
        else:
            nOfZList = [loadAscii(nOfZName, verbose=self.verbose) for nOfZName in nOfZNameList]
        
        zArr     = nOfZList[0][0]
        nOfZList = [nOfZ[1][:-1] for nOfZ in nOfZList]
        z_lower  = zArr[:-1]
        z_upper  = zArr[1:]
        z_middle = 0.5 * (z_lower + z_upper)
        kernel   = twopoint.NumberDensity(name, z_lower, z_middle, z_upper, nOfZList, ngal=nGalList, sigma_e=sigmaEpsList)
        return kernel
      
    def setKernels(self):
        labConv = LabelConvention()
        kernelN = self._makeKernel(labConv.lens, self.nOfZNameList[:self.nbTomoN], self.nGalList[:self.nbTomoN], None)
        kernelG = self._makeKernel(labConv.source, self.nOfZNameList[self.nbTomoN:], self.nGalList[self.nbTomoN:], self.sigmaEpsList)
        
        if kernelN is None:
            self.kernelList = [kernelG]
        elif kernelG is None:
            self.kernelList = [kernelN]
        else:
            self.kernelList = [kernelN, kernelG]
        return
    
    ## Mean
    def _loadNpyDataMat_Flinc(self, statsTag, randTag, verbose=True):
        statsList = statsTag.split('+')
        stock = []
        
        for stats in statsList:
            name = '%sdataMat_%s_%s_full.npy' % (self.prefix_Flinc, stats, randTag)
            data = np.load(name)
            stock.append(data)
            if verbose == True:
                print('Loaded \"%s\"' % name)
        
        ## Cut
        NList = [data.shape[0] for data in stock]
        N_min = min(NList)
        stock = [data[:N_min] for data in stock]
        stock = np.hstack(stock)
        
        if verbose == True:
            print('N = %s' % stock.shape[0])
        return stock
    
    def _makeMean_Flinc(self, statsTag, verbose=True):
        data = self._loadNpyDataMat_Flinc(statsTag, 'signal', verbose=verbose)
        mean = np.mean(data, axis=0)
        return mean
    
    def _interpolateCF(self, x, y):
        theta = self.theta
        inter = sp.interpolate.interp1d(np.log10(x), x*y, bounds_error=False, fill_value='extrapolate')
        CF    = inter(np.log10(theta)) / theta
        return theta, CF

    def makeTomoAngDict(self):
        labConv = LabelConvention()
        tomoAngDict = { ## Don't touch the order of this list
            labConv.w:       [self.__pairListNN, self.N_theta, self.theta],
            labConv.gamma_t: [self.__pairListNG, self.N_theta, self.theta], 
            labConv.gamma_x: [self.__pairListNG, self.N_theta, self.theta], 
            labConv.xi_p:    [self.__pairListGG, self.N_theta, self.theta], 
            labConv.xi_m:    [self.__pairListGG, self.N_theta, self.theta], 
            labConv.P_nn:    [self.__pairListNN, self.N_ell,   self.ell], 
            labConv.P_ne_E:  [self.__pairListNG, self.N_ell,   self.ell], 
            labConv.P_ne_B:  [self.__pairListNG, self.N_ell,   self.ell], 
            labConv.P_ee_E:  [self.__pairListGG, self.N_ell,   self.ell], 
            labConv.P_ee_B:  [self.__pairListGG, self.N_ell,   self.ell], 
            labConv.E_n:     [self.__pairListGG, self.nbModes, self.nArr], 
            labConv.B_n:     [self.__pairListGG, self.nbModes, self.nArr]
        }
        
        assert len(tomoAngDict) == len(labConv.kernelTypeDict)
        return tomoAngDict
    
    def _loadAsciiMean_CosmoSIS(self, stats, i, j, xArr, verbose=True):
        labConv = LabelConvention()
        
        if stats == labConv.w:
            name = '%sgalaxy_xi/bin_%d_%d.txt' % (self.prefix_CosmoSIS, j+1, i+1)
            yArr = loadAscii(name, verbose=verbose)
            theta, w = self._interpolateCF(xArr, yArr)
            return w
        
        if stats == labConv.gamma_t:
            name = '%sgalaxy_shear_xi/bin_%d_%d.txt' % (self.prefix_CosmoSIS, i+1, j+1)
            yArr = loadAscii(name, verbose=verbose)
            theta, gamma_t = self._interpolateCF(xArr, yArr)
            return gamma_t
          
        if stats == labConv.gamma_x:
            gamma_x = [0] * self.N_theta
            return gamma_x
        
        if stats == labConv.xi_p:
            name = '%spcfs/xi_p_bin_%d_%d.txt' % (self.prefix_CosmoSIS, j+1, i+1)
            xi_p = loadAscii(name, verbose=verbose)
            return xi_p
          
        if stats == labConv.xi_m:
            name = '%spcfs/xi_m_bin_%d_%d.txt' % (self.prefix_CosmoSIS, j+1, i+1)
            xi_m = loadAscii(name, verbose=verbose)
            return xi_m
          
        if stats == labConv.P_nn:
            name = '%sbandpower_galaxy/bin_%d_%d.txt' % (self.prefix_CosmoSIS, j+1, i+1)
            P_nn = loadAscii(name, verbose=verbose)
            return P_nn
        
        if stats == labConv.P_ne_E: #TODO
            P_ne_E = [0] * self.N_ell
            return P_ne_E
        
        if stats == labConv.P_ne_B:
            P_ne_B = [0] * self.N_ell
            return P_ne_B
        
        if stats == labConv.P_ee_E:
            name   = '%sbandpower_shear_e/bin_%d_%d.txt' % (self.prefix_CosmoSIS, j+1, i+1)
            P_ee_E = loadAscii(name, verbose=verbose)
            return P_ee_E
          
        if stats == labConv.P_ee_B:
            P_ee_B = [0] * self.N_ell
            return P_ee_B
          
        if stats == labConv.E_n:
            name = '%scosebis/bin_%d_%d.txt' % (self.prefix_CosmoSIS, j+1, i+1)
            E_n  = loadAscii(name, verbose=verbose)
            return E_n
        
        if stats == labConv.B_n:
            B_n = [0] * self.nbModes
            return B_n
        
        return None
    
    def _makeMean_CosmoSIS(self, statsTag, verbose=True):
        labConv = LabelConvention()
        tomoAngDict = self.makeTomoAngDict()
        statsList = statsTag.split('+')
        statsList_complete = tomoAngDict.keys()
        
        for stats in statsList:
            if stats not in statsList_complete:
                raise ValueError('\"%s\" not allowed' % statsTag)
        
        name  = '%sshear_xi_plus/theta.txt' % self.prefix_CosmoSIS
        xArr  = loadAscii(name, verbose=verbose) * (60.0 * 180.0 / np.pi) ## [arcmin]
        stock = []
        
        for stats, line in tomoAngDict.items():
            pairList = line[0]
          
            if stats in statsList:
                for i, j in pairList:
                    value = self._loadAsciiMean_CosmoSIS(stats, i, j, xArr, verbose=verbose)
                    stock.append(value)
        
        stock = np.concatenate(stock)
        return stock
    
    def makeMean(self, tag, name=None, statsTag=None, verbose=True):
        if tag == 'Flinc':
            return self._makeMean_Flinc(statsTag, verbose=verbose)
        if tag == 'CosmoSIS':
            return self._makeMean_CosmoSIS(statsTag, verbose=verbose)
        if tag == 'variable':
            return name ## consider "name" as the variable which contains the mean
        
        ## Otherwise, consider "name" as the path to the file
        ## which contains the mean data vector
        try:
            if name[-4:] == '.npy':
                mean = np.load(name).flatten()
            elif name[-4:] == '.fit' or name[-5:] == '.fits':
                mean = fits.getdata(name, 1).field(0)
            else:
                mean = np.loadtxt(name).flatten()
            if verbose:
                print('Loaded \"%s\"' % name)
        except:
            raise OSError('\"%s\" not found' % name)
        return mean
    
    ## Covariance
    def _makeCov_Flinc(self, statsTag, verbose=True):
        data = self._loadNpyDataMat_Flinc(statsTag, 'obs', verbose=verbose)
        cov  = np.cov(data, rowvar=0, ddof=1)
        return cov
    
    def _makePairIndex_2PCF(self, tomo1, tomo2, sign, order=-1):
        if order > 0:
            beginGG = 0
            beginNG = self.__nbPairsGG * 2
            beginNN = self.__nbPairsGG * 2 + self.__nbPairsNG
            split   = self.nbTomoG
            tomo1G  = tomo1 - 0
            tomo1N  = tomo1 - split
            tomo2G  = tomo2 - 0
            tomo2N  = tomo2 - split
            
            isGG    = asInt(tomo2 < split)
            isNG    = asInt((tomo1 < split) * (tomo2 >= split))
            isNN    = asInt(tomo1 >= split)
            ind     = isGG * (beginGG + pairIndex(self.nbTomoG, tomo1G, tomo2G) + self.__nbPairsGG * sign)
            ind    += isNG * (beginNG + tomo2N + self.nbTomoN * tomo1G) ## Vary lens tomo bins first
            ind    += isNN * (beginNN + pairIndex(self.nbTomoN, tomo1N, tomo2N))
          
        else:
            beginNN = 0
            beginNG = self.__nbPairsNN
            beginGG = self.__nbPairsNN + self.__nbPairsNG
            split   = self.nbTomoN
            tomo1N  = tomo1 - 0
            tomo1G  = tomo1 - split
            tomo2N  = tomo2 - 0
            tomo2G  = tomo2 - split
            
            isNN    = asInt(tomo2 < split)
            isNG    = asInt((tomo1 < split) * (tomo2 >= split))
            isGG    = asInt(tomo1 >= split)
          
        ind  = isNN * (beginNN + pairIndex(self.nbTomoN, tomo1N, tomo2N))
        ind += isNG * (beginNG + tomo2G + self.nbTomoG * tomo1N) ## Vary source tomo bins first
        ind += isGG * (beginGG + pairIndex(self.nbTomoG, tomo1G, tomo2G) + self.__nbPairsGG * sign)
        return ind
    
    def _makePairIndex_BP(self, tomo1, tomo2, order=-1):
        if order > 0:
            beginGG = 0
            beginNG = self.__nbPairsGG
            beginNN = self.__nbPairsGG + self.__nbPairsNG
            split   = self.nbTomoG
            tomo1G  = tomo1 - 0
            tomo1N  = tomo1 - split
            tomo2G  = tomo2 - 0
            tomo2N  = tomo2 - split
            
            isGG    = asInt(tomo2 < split)
            isNG    = asInt((tomo1 < split) * (tomo2 >= split))
            isNN    = asInt(tomo1 >= split)
            ind     = isGG * (beginGG + pairIndex(self.nbTomoG, tomo1G, tomo2G))
            ind    += isNG * (beginNG + tomo2N + self.nbTomoN * tomo1G) ## Vary lens tomo bins first
            ind    += isNN * (beginNN + pairIndex(self.nbTomoN, tomo1N, tomo2N))
          
        else:
            beginNN = 0
            beginNG = self.__nbPairsNN
            beginGG = self.__nbPairsNN + self.__nbPairsNG
            split   = self.nbTomoN
            tomo1N  = tomo1 - 0
            tomo1G  = tomo1 - split
            tomo2N  = tomo2 - 0
            tomo2G  = tomo2 - split
            
            isNN    = asInt(tomo2 < split)
            isNG    = asInt((tomo1 < split) * (tomo2 >= split))
            isGG    = asInt(tomo1 >= split)
            
        ind  = isNN * (beginNN + pairIndex(self.nbTomoN, tomo1N, tomo2N))
        ind += isNG * (beginNG + tomo2G + self.nbTomoG * tomo1N) ## Vary source tomo bins first
        ind += isGG * (beginGG + pairIndex(self.nbTomoG, tomo1G, tomo2G))
        return ind
    
    def _covListToMatrix_2PCF(self, data, cleanNaN=True, CTag='tot'):
        if cleanNaN:
            ind = np.isnan(data)
            data[ind] = 0.0
        
        tomoA1 = asInt(data[0]) - 1
        tomoA2 = asInt(data[1]) - 1
        tomoB1 = asInt(data[6]) - 1
        tomoB2 = asInt(data[7]) - 1
        signA  = asInt(data[4])
        signB  = asInt(data[10])
        binA   = asInt(data[5])
        binB   = asInt(data[11])
        if CTag == 'tot':
            value = data[12:].sum(axis=0)
        elif CTag == 'CNG2':
            value = data[12] + data[13] + data[14]
        elif CTag == 'CNG':
            value = data[12] + data[13]
        else:
            value = data[12]
        
        ## Make index
        indA  = self._makePairIndex_2PCF(tomoA1, tomoA2, signA, order=-1)
        indB  = self._makePairIndex_2PCF(tomoB1, tomoB2, signB, order=-1)
        indA  = binA + self.N_theta * indA
        indB  = binB + self.N_theta * indB
        d_tot = indA.max() + 1
        
        ## Fill the other triangle
        cov = np.zeros((d_tot, d_tot), dtype=float)
        cov[indA, indB] = value
        ind = np.arange(d_tot, dtype=int)
        cov[ind, ind] *= 0.5
        cov += cov.T
        return cov
    
    def _covListToMatrix_BP(self, data, cleanNaN=True, CTag='tot'):
        if cleanNaN:
            ind = np.isnan(data)
            data[ind] = 0.0
        
        tomoA1 = asInt(data[0]) - 1
        tomoA2 = asInt(data[1]) - 1
        tomoB1 = asInt(data[5]) - 1
        tomoB2 = asInt(data[6]) - 1
        binA   = asInt(data[4])
        binB   = asInt(data[9])
        value  = data[10]
        
        ## Make index
        indA  = self._makePairIndex_BP(tomoA1, tomoA2, order=-1)
        indB  = self._makePairIndex_BP(tomoB1, tomoB2, order=-1)
        indA  = binA + self.N_ell * indA
        indB  = binB + self.N_ell * indB
        d_tot = indA.max() + 1
        
        ## Fill the other triangle
        cov = np.zeros((d_tot, d_tot), dtype=float)
        cov[indA, indB] = value
        ind = np.arange(d_tot, dtype=int)
        cov[ind, ind] *= 0.5
        cov += cov.T
        return cov
    
    def _get2PCFIndexForCut(self, statsTag):
        labConv = LabelConvention()
        wTh = True if labConv.w in statsTag else False
        gT  = True if labConv.gamma_t  in statsTag else False
        xiP = True if labConv.xi_p in statsTag else False
        xiM = True if labConv.xi_m in statsTag else False
        ind = [wTh]*self.N_theta*self.__nbPairsNN + [gT]*self.N_theta*self.__nbPairsNG + [xiP]*self.N_theta*self.__nbPairsGG + [xiM]*self.N_theta*self.__nbPairsGG
        return ind
    
    def _getBPIndexForCut(self, statsTag):
        labConv = LabelConvention()
        statsList = statsTag.split('+')
        Pnn  = True if labConv.P_nn in statsTag else False
        PneE = True if labConv.P_ne_E in statsTag else False
        PneB = True if labConv.P_ne_B in statsTag else False
        PeeE = True if labConv.P_ee_E in statsTag else False
        PeeB = True if labConv.P_ee_B in statsTag else False
        
        if PneE or PeeE == True:
            ind = [Pnn]*self.N_ell*self.__nbPairsNN + [PneE]*self.N_ell*self.__nbPairsNG + [PeeE]*self.N_ell*self.__nbPairsGG
        else:
            ind = [Pnn]*self.N_ell*self.__nbPairsNN + [PneB]*self.N_ell*self.__nbPairsNG + [PeeB]*self.N_ell*self.__nbPairsGG
        return ind
    
    def _getCOSEBIIndexForCut(self, statsTag):
        labConv = LabelConvention()
        En  = True if labConv.E_n in statsTag else False
        Bn  = True if labConv.B_n in statsTag else False
        ind = [En]*self.nbModes*self.__nbPairsGG + [Bn]*self.nbModes*self.__nbPairsGG
        return ind
    
    @classmethod
    def getCategory(cls, statsTag):
        labConv = LabelConvention()
        statsList = statsTag.split('+')
        is2PCF    = False
        isBPE     = False
        isBPB     = False
        isCOSEBI  = False
        
        for stats in statsList:
            if stats in [labConv.w, labConv.gamma_t, labConv.gamma_x, labConv.xi_p, labConv.xi_m]:
                is2PCF   = is2PCF or True
            elif stats in [labConv.P_ne_E, labConv.P_ee_E]:
                isBPE    = isBPE or True
            elif stats in [labConv.P_ne_B, labConv.P_ee_B]:
                isBPB    = isBPB or True
            elif stats in [labConv.E_n, labConv.B_n]:
                isCOSEBI = isCOSEBI or True
        
        category = 1*int(is2PCF) + 2*int(isBPE) + 4*int(isBPB) + 8*int(isCOSEBI)
        return category

    def _makeCov_list(self, name, statsTag, cleanNaN=True, CTag='tot', verbose=True):
        covList  = loadAscii(name, verbose=verbose)
        category = TwoPointBuilder.getCategory(statsTag)
        
        if category == 1:
            cov = self._covListToMatrix_2PCF(covList, cleanNaN=cleanNaN, CTag=CTag)
            ind = self._get2PCFIndexForCut(statsTag)
        elif category in [0, 2, 4]:
            cov = self._covListToMatrix_BP(covList, cleanNaN=cleanNaN, CTag=CTag)
            ind = self._getBPIndexForCut(statsTag)
        elif category == 8:
            raise NotImplementedError('Reading COSEBI cov from list format not implemented')
        else:
            raise ValueError('statsTag = \"%s\" not allowed' % statsTag)
        
        cov = cov[ind].T[ind].T 
        return cov
    
    def makeCov(self, tag, name=None, statsTag=None, cleanNaN=True, CTag='tot', verbose=True):
        if tag == 'Flinc':
            return self._makeCov_Flinc(statsTag, verbose=verbose)
        if tag == 'list':
            return self._makeCov_list(name, statsTag, cleanNaN=cleanNaN, CTag=CTag, verbose=verbose)
        if tag == 'variable':
            return name ## consider "name" as the variable which contains the covariance
      
        ## Otherwise, consider "name" as the path to the file
        ## which contains the covariance matrix
        try:
            if name[-4:] == '.npy':
                cov = np.load(name)
            elif name[-4:] == '.fit' or name[-5:] == '.fits':
                cov = fits.getdata(name, 1)
            else:
                cov = np.loadtxt(name)
            if verbose:
                print('Loaded \"%s\"' % name)
        except:
            raise OSError('\"%s\" not found' % name)
        return cov
  
    ## Build up & save
    def makeTwoPoint_withCov(self, labConv, statsTag_c, mean, cov):
        tomoAngDict = self.makeTomoAngDict()
        statsList_c  = statsTag_c.split('+')
        statsList_c_complete = labConv.kernelTypeDict.keys()
      
        for stats_c in statsList_c:
            if stats_c not in statsList_c_complete:
                raise ValueError('\"%s\" not allowed' % statsTag_c)
    
        statsNameDict = {}
        builder = twopoint.SpectrumCovarianceBuilder()
        binInd  = 0
      
        for stats_c, line1, line2 in zip(labConv.kernelTypeDict.keys(), labConv.kernelTypeDict.values(), tomoAngDict.values()):
            ker1     = line1[0]
            ker2     = line1[1]
            type1    = line1[2]
            type2    = line1[3]
            #unit     = line1[4]
            pairList = line2[0]
            N_ang    = line2[1]
            angle    = line2[2]
            
            if stats_c in statsList_c:
                statsNameDict[(ker1, ker2, type1, type2)] = stats_c
                for i, j in pairList:
                    for angInd in range(N_ang):
                        x = angle[angInd]
                        y = mean[binInd]
                        binInd += 1
                        builder.add_data_point(ker1, ker2, type1, type2, i+1, j+1, x, angInd+1, y)
    
        ## Make TP
        builder.set_names(statsNameDict)
        spectra, cov_info = builder.generate(cov, 'arcmin')
        TP = TwoPointWrapper.from_spectra(spectra, kernels=self.kernelList, covmat_info=cov_info)
        #TP = TwoPointWrapper.from_spectra(spectra, kernels=self.kernelList, covmat_info=None)
        #TP = TwoPointWrapper.from_spectra(spectra, kernels=None, covmat_info=cov_info)
        #TP = TwoPointWrapper.from_spectra(spectra, kernels=None, covmat_info=None)
        return TP
  
    def makeTwoPoint_withoutCov(self, labConv, statsTag_c, mean):
        tomoAngDict = self.makeTomoAngDict()
        statsList_c  = statsTag_c.split('+')
        statsList_c_complete = labConv.kernelTypeDict.keys()
        
        for stats_c in statsList_c:
            if stats_c not in statsList_c_complete:
                raise ValueError('\"%s\" not allowed' % statsTag_c)
      
        spectra = []
        binInd  = 0
        
        for stats_c, line1, line2 in zip(labConv.kernelTypeDict.keys(), labConv.kernelTypeDict.values(), tomoAngDict.values()):
            ker1     = line1[0]
            ker2     = line1[1]
            type1    = line1[2]
            type2    = line1[3]
            unit     = line1[4]
            pairList = line2[0]
            N_ang    = line2[1]
            angle    = line2[2]
            
            if stats_c in statsList_c:
                sBuilder = SpectrumBuilder()
              
                for i, j in pairList:
                    value   = mean[binInd:binInd+N_ang]
                    binInd += N_ang
                    sBuilder.addTomo(i, j, angle, value)
          
                spec = sBuilder.makeSpectrum(stats_c, (type1, type2), unit, kernels=(ker1, ker2))
                #spec = sBuilder.makeSpectrum(stats_c, (type1, type2), unit)
                spectra.append(spec)
        
        ## Make
        TP = TwoPointWrapper.from_spectra(spectra, kernels=self.kernelList, covmat_info=None)
        #TP = TwoPointWrapper.from_spectra(spectra, kernels=None, covmat_info=None)
        return TP
  
###############################################################################
## Auxiliary functions

def loadAscii(name, sep=None, cmt='#', verbose=True):
    data = np.loadtxt(name, comments=cmt, delimiter=sep)
    if verbose is True:
        print('Loaded \"%s\"' % name)
    return data.T

def asInt(a, copy=True):
    if np.isscalar(a):
        return int(a)
    if type(a) is list:
        return np.array(a, dtype=int)
    return a.astype(int, copy=copy)
  
def pairIndex(N, i, j):
    ind  = (N + (N+1-i)) * i // 2
    ind += j - i
    return ind

###############################################################################
## Functions related to main functions

def defineScaleCuts():
    scDictDict = {
        'scale_cuts_none': {
            'use_stats': 'xiP xiM'
        },
      
        'scale_cuts_KV450': {
            'use_stats': 'xiP xiM',
            'keep_ang_xiP': '0.5 75',
            'keep_ang_xiM': '4 300'
        }
    }
    return scDictDict

def saveFitsTwoPoint(verbose=True):
    name = 'data/mockFootprint/otis/MFP_for_others/mean_sim_En+Bn_obs.dat'
    data1 = loadAscii(name)[:75]
    name = 'data/mockFootprint/otis/MFP_for_others/cov_th_En+Bn_obs.dat'
    data2 = loadAscii(name)[:75, :75]
    
    inputList = [
        ## statsTag, meanTag, meanName, covTag, covName
        ('gT+xiP+xiM', 'CosmoSIS', None, 'list', 'fakeCovariance.dat'),
        ('PneE+PeeE', 'Flinc', None, 'Flinc', None),
        ('En+Bn', 'file', 'data/mockFootprint/otis/MFP_for_others/mean_sim_En+Bn_obs.dat', 'file', 'data/mockFootprint/otis/MFP_for_others/cov_sim_En+Bn_obs.dat'),
        ('En', 'variable', data1, 'variable', data2),
        ('xiP+PeeE+En', 'Flinc', None, 'Flinc', None)
        
        #('PneE+PeeE', 'Flinc', None, 'Flinc', None)
    ]
    scOption = 'scale_cuts_none'
    
    ## Initialize instance
    labConv           = LabelConvention()
    scDictDict        = defineScaleCuts()
    statsList, scArgs = labConv.makeScaleCutsArgs(scDictDict[scOption])
    TPBuilder         = TwoPointBuilder()
    
    for statsTag_c, meanTag, meanName, covTag, covName in inputList:
        statsTag = labConv.customToDefaultStatsTag(statsTag_c)
        mean     = TPBuilder.makeMean(meanTag, name=meanName, statsTag=statsTag, verbose=verbose)
        cov      = TPBuilder.makeCov(covTag, name=covName, statsTag=statsTag, verbose=verbose)
        TP       = TPBuilder.makeTwoPoint_withCov(labConv, statsTag_c, mean, cov)
        #TP       = TPBuilder.makeTwoPoint_withoutCov(labConv, statsTag_c, mean)
        
        ## Scale cuts
        TP.cutScales(cutCross=scArgs[0], statsTag_tomoInd_tomoInd_list=scArgs[1], statsTag_binIndList_dict=scArgs[2])
        TP.keepScales(statsTag_tomoInd1_tomoInd2__angMin_angMax_dict=scArgs[3], statsTag__angMin_angMax_dict=scArgs[4])
        
        ## Save
        name = '%s_%s.fits' % (labConv.prefix, statsTag_c)
        TP.to_fits(name, overwrite=True, clobber=True)
        print('Saved \"%s\"' % name)
    return

def saveAsciiFakeCovariance():
    name = 'data/mockFootprint/milvus/cosmokids/thps_cov_kids1000_milvus_obs_complex_list.dat'
    data = loadAscii(name)
    
    ind = data[0] == 1
    data[0][~ind] += 1
    ind = data[1] == 1
    data[1][~ind] += 1
    ind = data[6] == 1
    data[6][~ind] += 1
    ind = data[7] == 1
    data[7][~ind] += 1
    
    name = 'fakeCovariance.dat'
    f = open(name, 'w')
    for line in data.T:
        f.write('%d  %d  %d  %d  %d  %d    %d  %d  %d  %d  %d  %d    % .10e  % .10e  % .10e\n' %\
            (line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14]))
    
    f.close()
    print('Saved \"%s\"' % name)
    return

def saveFitsTwoPoint_KV450(verbose=True):
    prefix = '../KCAP_scale_cuts/kv450_data/'
    name1  = '%sdata_vector/KV450_reweight_3x4x4_v2_good_xipm_mcor_5bin.dat' % prefix
    data1  = loadAscii(name1)
    data_p = data1[1:].T[:9].T.flatten()
    data_m = data1[1:].T[9:].T.flatten()
    mean   = np.concatenate([data_p, data_m])
    
    name2 = '%scov/cov_analytic.txt' % prefix
    cov   = loadAscii(name2)
    
    nOfZNameList = [
        '%snofz/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB01t03_blindA_Nz.asc' % prefix, 
        '%snofz/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB03t05_blindA_Nz.asc' % prefix, 
        '%snofz/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB05t07_blindA_Nz.asc' % prefix, 
        '%snofz/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB07t09_blindA_Nz.asc' % prefix, 
        '%snofz/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB09t12_blindA_Nz.asc' % prefix
    ]
    TPBuilder = TwoPointBuilder(
        nbTomoN=0, nbTomoG=5,
        nOfZNameList=nOfZNameList, 
        nGalList=[1]*5, 
        sigmaEpsList=[0.28]*5, 
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        prefix_Flinc=None, prefix_CosmoSIS=None
    )
    scOption = 'scale_cuts_none'
    
    ## Initialize instance
    labConv           = LabelConvention()
    scDictDict        = defineScaleCuts()
    statsList, scArgs = labConv.makeScaleCutsArgs(scDictDict[scOption])
    TP = TPBuilder.makeTwoPoint_withCov(labConv, '+'.join(statsList), mean, cov)
    
    ## Scale cuts
    TP.cutScales(cutCross=scArgs[0], statsTag_tomoInd_tomoInd_list=scArgs[1], statsTag_binIndList_dict=scArgs[2])
    TP.keepScales(statsTag_tomoInd1_tomoInd2__angMin_angMax_dict=scArgs[3], statsTag__angMin_angMax_dict=scArgs[4])
    
    ## Save
    name = '%s_%s.fits' % (labConv.prefix, 'KV450')
    TP.to_fits(name, overwrite=True, clobber=True)
    if verbose:
        print('Saved \"%s\"' % name)
    return

def saveFitsTwoPoint_fakeCTerm(verbose=True):
    prefix = '../KCAP_scale_cuts/kv450_data/'
    name3  = '%ssystematics/KV450_ALL_c12_treecorr.out' % prefix
    data3  = loadAscii(name3)
    cTerm  = data3[3]
    
    tpType4 = twopoint.Types.galaxy_position_fourier    ## GPF
    tpType5 = twopoint.Types.galaxy_shear_emode_fourier ## GEF
    tpType6 = twopoint.Types.galaxy_shear_bmode_fourier ## GBF
    TPBuilder = TwoPointBuilder(verbose=verbose)
    
    sBuilder = SpectrumBuilder()
    for i in range(5):
        for j in range(i, 5):
            sBuilder.addTomo(i, j, TPBuilder.theta, cTerm)
    spec1 = sBuilder.makeSpectrum('2d_cterm', (tpType4, tpType4), None)
    
    sBuilder = SpectrumBuilder()
    for i in range(5):
        for j in range(i, 5):
            sBuilder.addTomo(i, j, TPBuilder.theta, cTerm)
    spec2 = sBuilder.makeSpectrum('cterm_cos', (tpType5, tpType5), None)
    
    sBuilder = SpectrumBuilder()
    for i in range(5):
        for j in range(i, 5):
            sBuilder.addTomo(i, j, TPBuilder.theta, cTerm)
    spec3 = sBuilder.makeSpectrum('cterm_sin', (tpType6, tpType6), None)
    
    TP = TwoPointWrapper.from_spectra([spec1, spec2, spec3], kernels=None, covmat_info=None)
    
    name = 'twoPoint_fakeCTerm.fits'
    TP.to_fits(name, overwrite=True, clobber=True)
    if verbose:
        print('Saved \"%s\"' % name)
    return

###############################################################################
## Check

def initialize():
    labConv = LabelConvention()
    name = '%s_PneE+PeeE.fits' % labConv.prefix
    return name

def printTwoPointHDU(ind=1):
    name = initialize()
    hdr  = fits.getheader(name, ind)
    data = fits.getdata(name, ind)
    
    print()
    print(hdr.tostring(sep='\n'))
    print(data)
    return 

def printTwoPoint_fromFile():
    name    = initialize()
    HDUList = fits.open(name)
    
    ## Check default HDU
    hdr0    = HDUList[0].header
    print()
    print('Check default HDU:')
    
    try:
        dummy   = hdr0['SIMPLE']
        HDUList = HDUList[1:]
        nbHDU   = len(HDUList)
        print('  Passed.')
    except:
        print('  Failed but continue.')
        print('  Keep in mind it means that this file was not generated in')
        print('  the standard way.')
    
    hdrList  = [HDU.header for HDU in HDUList]
    dataList = [HDU.data for HDU in HDUList]
    
    print()
    checkCovariance_fromFile(hdrList)
    print()
    checkSpectra_fromFile(hdrList)
    print()
    checkKernels_fromFile(hdrList, dataList)
    print()
    checkNGal_fromFile(hdrList)
    
    #print()
    #print()
    #print(hdrList[HDUInd-1].tostring(sep='\n'))
    return 

def checkCovariance_fromFile(hdrList):
    try:
        extNameList = [hdr['EXTNAME'] for hdr in hdrList]
    except:
        print('  One of the extensions doesn\'t have \"EXTNAME\" key.')
        print('  That is not normal. Exit.')
        return
    
    print('Check covariance:')
    
    try:
        ind = extNameList.index('COVMAT')
        print('  Passed.')
        hdr   = hdrList[ind]
        count = 0
        stock = []
        print()
        print('Dimension of the covariance:')
        
        while count >= 0:
            try:
                stock.append(hdr['STRT_%d' % count])
                count += 1
            except:
                count = -1
        count = len(stock)
        stock.append(hdr['NAXIS1'])
        
        for i in range(count):
            print('  %8s = %3d' % (hdr['NAME_%d' % i], stock[i+1]-stock[i]))
        print('     Total = %3d' % stock[-1])

    except:
        print('  Failed. Skipped.')
    return

def checkSpectra_fromFile(hdrList):
    print('Dimension of the data vector:')
    print('     stats   dim  tomo_1  tomo_2  N_ang')
    
    for hdr in hdrList:
        try:
          print('  %8s   %3d      %2d      %2d     %2d' % (hdr['EXTNAME'], hdr['NAXIS2'], hdr['N_ZBIN_1'], hdr['N_ZBIN_2'], hdr['N_ANG']))
        except:
          continue
    return

def checkKernels_fromFile(hdrList, dataList):
    print('Redshift ranges:')
    print('     kernel     z_min     z_max')
    
    for hdr, data in zip(hdrList, dataList):
        try:
            print('  %9s  %8f  %8f' % (hdr['EXTNAME'], data['Z_LOW'][0], data['Z_HIGH'][-1]))
        except:
            continue
    return 

def checkNGal_fromFile(hdrList):
    print('Galaxy number densities:')
    
    for hdr in hdrList:
        try:
            dummy = hdr['NGAL_1']
            count = 1
            print('  %s' % hdr['EXTNAME'])
            while count > 0:
                try:
                    print('    NGAL%d = %8f' % (count, hdr['NGAL_%d' % count]))
                    count += 1
                except:
                    count = -1
        except:
            continue
    return

def loadFitsTwoPoint():
    name = initialize()
    try:
        TP = twopoint.TwoPointFile.from_fits(name, covmat_name='COVMAT')
    except:
        TP = twopoint.TwoPointFile.from_fits(name, covmat_name=None)
    printTwoPoint_fromObj(TP)
    return

def printTwoPoint_fromObj(TP):
    print('Spectra:')
    print()
    for spectrum in TP.spectra:
        printSpectrum_fromObj(spectrum)
        print()
    
    print('Kernels:')
    print()
    for kernel in TP.kernels:
        printKernel_fromObj(kernel)
        print()
    
    ##print(TP.windows)
    print('Covariance:')
    print()
    if hasattr(TP, 'covmat_info') and TP.covmat_info is not None:
        printCovMatInfo_fromObj(TP.covmat_info)
    ##print(TP._spectrum_index)
    return

def printSpectrum_fromObj(spectrum):
    print('name = %s' % spectrum.name)
    print('tomoInd1 = %s' % spectrum.bin1)
    print('tomoInd2 = %s' % spectrum.bin2)
    print('pair = %s' % spectrum.bin_pairs)
    ##print(spectrum.type1)
    ##print(spectrum.type2)
    print('nOfZ1 = %s' % spectrum.kernel1)
    print('nOfZ2 = %s' % spectrum.kernel2)
    print('angInd = %s' % spectrum.angular_bin)
    #print('ang = %s' % spectrum.angle)
    ##print(spectrum.angle_min)
    ##print(spectrum.angle_max)
    #print('mean = %s' % spectrum.value)
    ##print(spectrum.npairs)
    ##print(spectrum.varxi)
    ##print(spectrum.windows)
    #print('std = %s' % spectrum.error)
    ##print(spectrum.metadata)
    #print('unit = %s' % spectrum.angle_unit)
    print('binInd = %s' % spectrum.dv_index)
    ##print(spectrum.extra_cols)
    return

def printKernel_fromObj(kernel):
    print('name = %s' % kernel.name)
    print('zlow = %s' % kernel.zlow)
    print('z = %s' % kernel.z)
    print('zhigh = %s' % kernel.zhigh)
    print('nbNOfZ = %s' % kernel.nbin)
    print('nbZBins = %s' % kernel.nsample)
    print('nArr = %s' % kernel.nzs)
    print('n_gal = %s' % kernel.ngal)
    print('sigma_e = %s' % kernel.sigma_e)
    return

def printCovMatInfo_fromObj(covMatInfo):
    pass
    #print('name = %s' % covMatInfo.name)
    print('names = %s' % covMatInfo.names)
    print('lengths = %s' % covMatInfo.lengths)
    print('starts = %s' % covMatInfo.starts)
    print('covmat = %s' % covMatInfo.covmat)
    #print('diagonal = %s' % covMatInfo.diagonal)
    return
  
###############################################################################

