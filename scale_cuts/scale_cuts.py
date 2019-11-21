from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.datablock.cosmosis_py import errors
import numpy as np
#import twopoint
import wrapper_twopoint as wtp


def deleteComment(line):
    ind = line.find(';')
    if ind < 0:
        return line
    return line[:ind]

def loadScaleCutsFile(name):
    f = open(name, 'r')
    stock = [line for line in f]
    f.close()
    
    indList = [i for i, line in enumerate(stock) if line[0] == '[']
    indList.append(len(stock))
    indList1 = indList[:-1]
    indList2 = indList[1:]
    scDictDict = {}
    
    for ind1, ind2 in zip(indList1, indList2):
        key   = stock[ind1].strip(' []\n')
        block = stock[ind1+1:ind2]
        block = [deleteComment(line) for line in block] ## Delete comments
        block = [line.split('=') for line in block]
        block = [[word.strip(' \n') for word in line] for line in block]
        block = [line for line in block if len(line) > 1]
        scDict = {}
        for line in block:
            scDict[line[0]] = line[1]
        scDictDict[key] = scDict
    return scDictDict

def setup(options):
    config = {}
    
    ## Read inputs from datablock
    config['output_section_name'] = options.get_string(option_section, 'output_section_name', default='likelihood')
    try:
        config['data_and_covariance_fits_filename'] = options.get_string(option_section, 'data_and_covariance_fits_filename')
    except errors.BlockNameNotFound:
        raise NameError('data_and_covariance_fits_filename cannot be empty')
    
    ## Read extension names for data outputs
    config['xi_plus_extension_name']                  = options.get_string(option_section, 'xi_plus_extension_name', default='xiP')
    config['xi_minus_extension_name']                 = options.get_string(option_section, 'xi_minus_extension_name', default='xiM')
    config['bandpower_ggl_extension_name']            = options.get_string(option_section, 'bandpower_ggl_extension_name', default='PneE')
    config['bandpower_e_cosmic_shear_extension_name'] = options.get_string(option_section, 'bandpower_e_cosmic_shear_extension_name', default='PeeE')
    config['cosebis_extension_name']                  = options.get_string(option_section, 'cosebis_extension_name', default='En')
    
    ## Read section names for theory outputs
    config['xi_plus_section_name']                  = options.get_string(option_section, 'xi_plus_section_name', default='shear_xi_plus')
    config['xi_minus_section_name']                 = options.get_string(option_section, 'xi_minus_section_name', default='shear_xi_minus')
    config['bandpower_ggl_section_name']            = options.get_string(option_section, 'bandpower_ggl_section_name', default='bandpower_ggl')
    config['bandpower_e_cosmic_shear_section_name'] = options.get_string(option_section, 'bandpower_e_cosmic_shear_section_name', default='bandpower_e_cosmic_shear')
    config['cosebis_section_name']                  = options.get_string(option_section, 'cosebis_section_name', default='cosebis')
    
    ## Read scale cuts
    try:
        config['scale_cuts_filename'] = options.get_string(option_section, 'scale_cuts_filename')
    except errors.BlockNameNotFound:
        raise NameError('scale_cuts_filename cannot be empty')
    config['scale_cuts_option'] = options.get_string(option_section, 'scale_cuts_option', default='scale_cuts_fiducial')
    
    ## Load data & cov file
    try:
        TP_data = wtp.TwoPointWrapper.from_fits(config['data_and_covariance_fits_filename'], covmat_name='COVMAT')
    except:
        raise OSError('\"%s\" not found' % config['data_and_covariance_fits_filename'])
    
    ## Load scale cuts file
    try:
        scDictDict = loadScaleCutsFile(config['scale_cuts_filename'])
    except:
        raise OSError('\"%s\" not found' % config['scale_cuts_filename'])
    
    ## Make scale cuts arguments
    try:
        scDict = scDictDict[config['scale_cuts_option']]
    except:
        raise KeyError('Bad scale cuts option: \"%s\"' % config['scale_cuts_option'])
    labConv = wtp.LabelConvention(xi_p=config['xi_plus_extension_name'],
                                  xi_m=config['xi_minus_extension_name'],
                                  P_ne_E=config['bandpower_ggl_extension_name'],
                                  P_ee_E=config['bandpower_e_cosmic_shear_extension_name'],
                                  E_n=config['cosebis_extension_name'])
    statsList, scArgs = labConv.makeScaleCutsArgs(scDict) ## Here, we convert the keys from the default ones to the custom ones.
    config['scale_cuts_arguments'] = scArgs
    config['label_convention']     = labConv
    
    ## Do scale cuts to data and cov
    TP_data.cutScales(cutCross=scArgs[0], statsTag_tomoInd_tomoInd_list=scArgs[1], statsTag_binIndList_dict=scArgs[2])
    TP_data.keepScales(statsTag_tomoInd1_tomoInd2__angMin_angMax_dict=scArgs[3], statsTag__angMin_angMax_dict=scArgs[4])
    print('  Did scale cuts to data & cov')
    
    ## Don't put this before scale cuts, because TP_data.choose_data_sets does 
    ## not modify the covmat_info attribute but only the covmat attribute.
    ## While one can cut a correctly cut covariance, the scale cuts on this
    ## covariance, which relies on covmat_info, can result in errors.
    statsList_c           = [labConv.defaultToCustomStatsTag(stats) for stats in statsList]
    config['use_stats']   = statsList
    config['use_stats_c'] = statsList_c
    TP_data.choose_data_sets(statsList_c)
    
    ## Extract the vector & matrix & put in config dict
    config['data']       = TP_data.makeMeanVector()
    config['covariance'] = TP_data.covmat
    
    config['simulate'] = options.get_bool(option_section, 'simulate', default=False)
    config['mock_filename'] = options.get_string(option_section, 'mock_filename', default="")
    ## Test
    #wtp.printTwoPoint_fromObj(TP_data)
    return config

def execute(block, config):
    print('Gathering theory outputs to make a vector')
    
    ## Save data and cov in the data block
    output_section_name = config['output_section_name']
    block[output_section_name, 'data']       = config['data']
    block[output_section_name, 'covariance'] = config['covariance']
    
    ## Define some things first
    labConv    = config['label_convention']
    nbTomo_max = 1000
    spectra    = []
    
    ## Don't change the order of this list
    ## Read as: [section_name, extension_name, angle_name, isGGL]
    sectionNameList = [
        [config['xi_plus_section_name'],                  config['xi_plus_extension_name'],                  'theta_bin_1_1', False],
        [config['xi_minus_section_name'],                 config['xi_minus_extension_name'],                 'theta_bin_1_1', False],
        [config['bandpower_ggl_section_name'],            config['bandpower_ggl_extension_name'],            'ell',           True],
        [config['bandpower_e_cosmic_shear_section_name'], config['bandpower_e_cosmic_shear_extension_name'], 'ell',           False],
        [config['cosebis_section_name'],                  config['cosebis_extension_name'],                  'cosebis_n',     False],
    ]
    
    for line in sectionNameList:
        section_name   = line[0]
        extension_name = line[1]
        angle_name     = line[2]
        isGGL          = line[3]
      
        if extension_name not in config['use_stats']:
            print('  Skipped %s' % extension_name)
            continue
        
        elif not block.has_section(section_name):
            raise AssertionError('You specified \"%s\" in use_stats. When I was computing for \"%s\", I tried to\n\
                grab theory values from \"%s\" but could not find anything.' % (config['use_stats'], extension_name, section_name))
        
        ## Initialize!
        sBuilder = wtp.SpectrumBuilder()
        
        ## Is it BP? If yes, then ell needs to be calculated instead of calling directly from data block.
        if angle_name == 'ell':
            lower = block[section_name, 'l_min_vec']
            upper = block[section_name, 'l_max_vec']
            angle = 0.5 * (np.log10(lower) + np.log10(upper))
            angle = 10**angle
        else:
            angle = block[section_name, angle_name]
        
        ## Is it GGL? If yes, the indices i & j run differently from the others.
        if isGGL:
            for i in range(nbTomo_max):
                if not block.has_value(section_name, 'bin_%d_%d' % (i+1, 1)):
                    print('  %s stops at lens bin %d' % (section_name, i))
                    break
          
                for j in range(nbTomo_max):
                    value_name  = 'bin_%d_%d' % (i+1, j+1)
                    try:
                        value = block[section_name, value_name]
                        sBuilder.addTomo(i, j, angle, value)
                    except:
                        break
        else:
            for i in range(nbTomo_max):
                if not block.has_value(section_name, 'bin_%d_%d' % (i+1, i+1)):
                    print('  %s stops at bin %d' % (section_name, i))
                    break
                
                for j in range(i, nbTomo_max):
                    value_name  = 'bin_%d_%d' % (j+1, i+1)
                    try:
                        value = block[section_name, value_name]
                        sBuilder.addTomo(i, j, angle, value)
                    except:
                        break
      
        ## These are some necessary components for generating the mean vector (spectrum).
        line  = labConv.kernelTypeDict[extension_name]
        ker1  = line[0]
        ker2  = line[1]
        type1 = line[2]
        type2 = line[3]
        unit  = line[4]
        
        spec  = sBuilder.makeSpectrum(extension_name, (type1, type2), unit, kernels=(ker1, ker2))
        spectra.append(spec)
      
    ## Make. Ta da!
    TP_theory = wtp.TwoPointWrapper.from_spectra(spectra, kernels=None, covmat_info=None)
    #TP_theory = wtp.TwoPointWrapper.from_fits('twoPoint_PneE+PeeE.fits', covmat_name=None) ## For test
    
    ## Do scale cuts to theory
    scArgs = config['scale_cuts_arguments']
    TP_theory.cutScales(cutCross=scArgs[0], statsTag_tomoInd_tomoInd_list=scArgs[1], statsTag_binIndList_dict=scArgs[2])
    TP_theory.keepScales(statsTag_tomoInd1_tomoInd2__angMin_angMax_dict=scArgs[3], statsTag__angMin_angMax_dict=scArgs[4])
    print('  Did scale cuts to theory')
    
    ## Extract the vector & put in block as output_section_name
    block[output_section_name, 'theory'] = TP_theory.makeMeanVector()

    if config["simulate"]:
        mu = block[output_section_name, 'theory']
        cov = block[output_section_name, 'covariance']
        s = np.random.multivariate_normal(mu, cov)
        TP_theory.replaceMeanVector(s)
        if config["mock_filename"] != "":
            TP_theory.to_fits(config["mock_filename"])
    
    ## Some print functions for debug
    #print()
    #print(block.keys())
    #print(block.sections())
    #print(block['nz_kv450_5bin', 'nbin'])
    #print()
    
    #wtp.printTwoPoint_fromObj(TP_theory)
    #print(block[output_section_name, 'data'].shape)
    #print(block[output_section_name, 'theory'].shape)
    #print(block[output_section_name, 'covariance'].shape)
    return 0
