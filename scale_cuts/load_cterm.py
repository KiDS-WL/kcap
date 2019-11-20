from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.datablock.cosmosis_py import errors
import numpy as np
#import twopoint
import wrapper_twopoint as wtp


def setup(options):
    config = {}
    
    ## Read inputs from datablock
    try:
        config['cterm_filename'] = options.get_string(option_section, 'cterm_filename')
        print('cterm_filename = %s' % config['cterm_filename'])
        
        config['cterm_extention_names']    = options.get_string(option_section, 'cterm_extention_names', default='2d_cterm cterm_cos cterm_sin')
        print('cterm_extention_names = %s' % config['cterm_extention_names'])
        cterm_extention_name_list          = [str_.strip() for str_ in config['cterm_extention_names'].split()]
        config['2d_cterm_extention_name']  = cterm_extention_name_list[0]
        config['cterm_cos_extention_name'] = cterm_extention_name_list[1]
        config['cterm_sin_extention_name'] = cterm_extention_name_list[2]
        
        #output_folder_names = options.get_string(option_section, 'output_folder_names', default='2d_cterm cterm_cos cterm_sin')
        #output_folder_names = [p.strip().lower() for p in output_folder_names.split()]
        #print('output_folder_names = %s' % output_folder_names)

    except errors.BlockNameNotFound:
        print("***No cterm file was given***\n")
        return config

    ## Load c-term file
    try:
        TP_cTerm = wtp.TwoPointWrapper.from_fits(config['cterm_filename'], covmat_name=None)
    except:
        raise OSError('\"%s\" not found' % config['cterm_filename'])
    
    ## Put in config
    spec_2d  = TP_cTerm.get_spectrum(config['2d_cterm_extention_name'])
    spec_cos = TP_cTerm.get_spectrum(config['cterm_cos_extention_name'])
    spec_sin = TP_cTerm.get_spectrum(config['cterm_sin_extention_name'])
    config['unique_pairs'] = spec_2d.get_bin_pairs()
    
    for i, j in config['unique_pairs']:
        theta, value_2d  = spec_2d.get_pair(i, j)
        theta, value_cos = spec_cos.get_pair(i, j)
        theta, value_sin = spec_sin.get_pair(i, j)
        config['2d_cterm_%d_%d' % (i, j)] = value_2d
        config['cterm_cos_%d_%d' % (i, j)] = value_cos
        config['cterm_sin_%d_%d' % (i, j)] = value_sin
    
    config['theta'] = theta
    
    ## Save here?
    
    #print()
    #print(block.keys())
    #print(block.sections())
    #print(block['nz_kv450_5bin', 'nbin'])
    #print()
    
    return config


def execute(block, config):
    
    ## No c-term file
    if len(config) == 0:
        return 0
  
    ## Save everything in the data block
    block[config['2d_cterm_extention_name'],  'theta'] = config['theta']
    block[config['cterm_cos_extention_name'], 'theta'] = config['theta']
    block[config['cterm_sin_extention_name'], 'theta'] = config['theta']
    
    for i, j in config['unique_pairs']:
        value_name = 'bin_%d_%d' % (j, i)
        block[config['2d_cterm_extention_name'],  value_name] = config['2d_cterm_%d_%d' % (i, j)]
        block[config['cterm_cos_extention_name'], value_name] = config['cterm_cos_%d_%d' % (i, j)]
        block[config['cterm_sin_extention_name'], value_name] = config['cterm_sin_%d_%d' % (i, j)]

    return 0
