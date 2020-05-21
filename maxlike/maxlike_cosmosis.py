from __future__ import print_function
import numpy as np
from cosmosis.runtime.pipeline import LikelihoodPipeline
from cosmosis.runtime.config import Inifile
import numpy as np
from math import pi
from scipy.optimize import minimize
import os
from chain_functions import find_cols_for_params, input_names, parameter_names
from argparse import ArgumentParser


# for testing this code
# foldername="/Users/marika_asgary/Documents/CosmicShear/repos/kids1000_chains/first_chains/impact_of_nuisance_parameters/free_m_bias_correlated/cosebis/"
# valuefile = foldername+"/config/values.ini"
# inifile = foldername+"/config/pipeline.ini"
# multinest_file = foldername+"/chain/output_multinest_A.txt"

# print('using the multinest file:'+multinest_file+' to set the starting point for maxlike')
# file=open(multinest_file)
# line=file.readline()
# parameter_list=(line.replace('#','')).split()
# multi_chain = np.loadtxt(multinest_file,comments='#')
# cols, param_names=find_cols_for_params(parameter_list,input_names,parameter_names,['post','like'])
# max_ind = np.argmax(multi_chain[:,cols[0]])
# best_fit = multi_chain[max_ind,:]


parser = ArgumentParser(description='runs maxlike chains using an input ini file')
parser.add_argument("-i", "--inifile", dest="inifile",
    help="Full Input ini file name and path", metavar="inifile",required=True)

parser.add_argument("-m", "--multinest_start_point", dest="multinest_file",
    help="multinest file from which the starting point will be taken", metavar="multinest_file",required=False)

parser.add_argument("-o", "--output_name", dest="output_name",
    help="name of the output file", metavar="output_name",required=True)

parser.add_argument("--max_post",action="store_true",
    help="Do you want max_post? If not set uses max_like",required=False)

parser.add_argument("-s","--ouput_steps",action="store_true",
    help="Output all steps of finding the minimum, default is false",required=False)

parser.add_argument("--maxiter", dest='maxiter',default = 1600, type=int,
    help="maximum number of iterations",required=False)

parser.add_argument("--best_fit_value", dest="best_fit_value",
    help="If set makes a values file with the best fit from the minimiser and the priors taken from the input values file, the name of the file needs to be given as a string.", metavar="best_fit_value",required=False)

parser.add_argument("--best_fit_priors", dest="best_fit_priors",
    help="If set makes a priors file with the best fit from the minimiser and the priors taken from the input priors file, the name of the file needs to be given as a string.", metavar="best_fit_priors. If best_fit_value given it will write this file but call it what ever name you gave to prior_best_fit_value",required=False)

args        = parser.parse_args()
inifile     = args.inifile
output_name = args.output_name
maxiter     = args.maxiter


if(args.best_fit_value):
  if(args.best_fit_priors):
    output_prior_file_name = args.best_fit_priors
  else:
    output_prior_file_name = "prior_"+args.best_fit_value


ini = Inifile(inifile) 
extra_params_ini=ini.get("pipeline","extra_output").split()

extra_params=[]
for param in extra_params_ini:
  param=param.replace("/","--")
  extra_params.append(param)

# run setup
pipeline = LikelihoodPipeline(ini) 
# list of varied parameter names
param_names_varied = pipeline.varied_params

if args.multinest_file:
  # get the starting point from the best fit of the multinest chain if the file is given
  multinest_file=args.multinest_file
  print('using the multinest file:'+multinest_file+' to set the starting point for maxlike')
  file=open(multinest_file)
  line=file.readline()
  parameter_list=(line.replace('#','')).split()
  multi_chain = np.loadtxt(multinest_file,comments='#')
  cols, param_names=find_cols_for_params(parameter_list,input_names,parameter_names,['post','like'])
  if(args.max_post):
    max_ind = np.argmax(multi_chain[:,cols[0]])
  else:
    max_ind = np.argmax(multi_chain[:,cols[1]])
  best_fit = multi_chain[max_ind,:]
  params_fiducial = best_fit[0:len(param_names_varied)]
  str_start = 'best fit of multinest file: '+multinest_file
  start_point_str = 'start values:'
  for param in best_fit:
    start_point_str+=' '+'%.5f' % param
else:
  # Otherwise use the fiducial values from the values.ini file:
  params_fiducial = pipeline.start_vector()
  str_start = 'fiducial start taken from the values file: '+ini.get("pipeline","values")
  start_point_str = 'start values:'
  for param in best_fit:
    start_point_str+=' '+'%.5f' % param



# set the pipeline properties
pipeline.quiet = True
pipeline.debug = False
pipeline.timing = False

# the normalised vector is used by the minimiser
start_vector  = pipeline.normalize_vector(params_fiducial)
params_values = pipeline.denormalize_vector(start_vector)

def run(normalised_values):
  if((normalised_values>1.).any() or (normalised_values<0.).any()):
    return np.inf
  else:
    params_values = pipeline.denormalize_vector(normalised_values)
    post, extra   = pipeline.posterior(params_values)
    prior         = pipeline.prior(params_values)
    like = post - prior
    if(args.ouput_steps):
      for p in params_values:
        if(p>0.01):
          file.write('%.5f' % p)
        else:
          file.write('%.5e' % p)
        file.write(" ")
      for p in extra:
        if(p>0.01):
          file.write('%.5f' % p)
        else:
          file.write('%.5e' % p)
        file.write(" ")
      # 
      file.write('%.5f' % like)
      file.write(" ")
      file.write('%.5f' % post)
      file.write("\n")
    if (args.max_post):
      return -2.*post
    else:
      return -2.*like

min_bound = np.asarray([0.0 for p in pipeline.varied_params])
max_bound = np.asarray([1.0 for p in pipeline.varied_params])
# have to normalise the parameters to go between -1 and 1 or something
nparam = len(pipeline.varied_params)


##############################################################################################

print("\n\n\n nelder mead \n\n")

file = open(output_name, 'w')
file.write("#")
for param in pipeline.varied_params:
  file.write(str(param))
  file.write(" ")

for param in extra_params:
  print(param)
  file.write(str(param))
  file.write(" ")

file.write("like post\n")
file.write("## maxlike sampler\n")
str_max_post = "Using maximum "
if(args.max_post):
  str_max_post+= " posterior values"
else:
  str_max_post+= " likelihood values"

file.write("## "+str_max_post+"\n")
file.write("## maxiter = "+str(maxiter)+"\n")
file.write("## input taken from "+str_start+"\n")
file.write("## "+start_point_str+"\n")
file.write("## START_OF_PARAMS_INI\n")

# now write in the ini file configuration
with open(inifile) as fp:
  line = fp.readline()
  while line:
    if(line[0]!=";"):
      # print('#'+line.strip())
      file.write("## ")
      file.write(line)
    line = fp.readline()

file.write("## END_OF_PARAMS_INI\n")
file.write("## START_OF_VALUES_INI\n")
valuefile = ini.get("pipeline","values")
with open(valuefile) as fp:
  line = fp.readline()
  while line:
    if(line[0]!=";"):
      # print('#'+line.strip())
      file.write("## ")
      file.write(line)
    line = fp.readline()

file.write("## END_OF_VALUES_INI\n")

try:
  priorfile = ini.get("pipeline","priors")
  file.write("## START_OF_PRIORS_INI\n")
  with open(priorfile) as fp:
    line = fp.readline()
    while line:
      if(line[0]!=";"):
        # print('#'+line.strip())
        file.write("## ")
        file.write(line)
      line = fp.readline()
  file.write("## END_OF_PRIORS_INI\n")
except:
  print("no prior file given")



res_nelder_mead = minimize(run, start_vector, method='Nelder-Mead', 
          options={'maxiter':maxiter,'xatol':1e-3, 'fatol':1e-2,'adaptive':True, 'disp':True})

output_params = pipeline.denormalize_vector(res_nelder_mead.x)
post, extra   = pipeline.posterior(output_params)
prior         = pipeline.prior(output_params)
like          = post - prior
for p in output_params:
  if(p>0.01):
    file.write('%.5f' % p)
  else:
    file.write('%.5e' % p)
  file.write(" ")

for p in extra:
  if(p>0.01):
    file.write('%.5f' % p)
  else:
    file.write('%.5e' % p)
  file.write(" ")

file.write('%.5f' % like)
file.write(" ")
file.write('%.5f' % post)
file.write("\n")

file.write("# post_min=")
file.write(str(res_nelder_mead.fun/2.))
file.write("# n_eval=")
file.write(str(res_nelder_mead.nfev))
file.write("\n")
if res_nelder_mead.success:
  file.write("# Success\n")
else:
  file.write("# Failure:"+str(res_nelder_mead.status)+"\n")

file.write("# ")
file.write(res_nelder_mead.message)
file.write("\n")

file.close()

##############################################################################################

print("Nelder-Mead results:")
output_params=pipeline.denormalize_vector(res_nelder_mead.x)
p=0
for param in pipeline.varied_params:
  print(param)
  print(param,output_params[p],params_fiducial[p])
  p+=1

print("min chi2=",res_nelder_mead.fun)


# now making an output values file 
if args.best_fit_value:
  filename = args.best_fit_value
  pipeline.create_ini(output_params,filename)
# if there is a prior file then change the values for the mean there as well.
  try:
    priorfile = ini.get("pipeline","priors")
    import collections
    output = collections.defaultdict(list)
    with open(priorfile) as fp:
      line = fp.readline()
      while line:
        if(line[0]=='['):
          section_name = line.strip()[1:-1]
          print(section_name)
          # what prior type is used? 
        elif(line.split()):
          if((line.split()[0]!=";")):
            param_name = line.split()[0]
            prior_type = line.split()[2]
            if(prior_type=='gaussian'):
              sigma      = float(line.split()[4])
              for param, x in zip(pipeline.varied_params, output_params):
                if((param.section == section_name) & (param.name == param_name)):
                  output[section_name].append("%s  =  %s    %r    %r\n"% (param.name,prior_type,x,sigma))
        line = fp.readline()
    file = open(output_prior_file_name, 'w')
    for section, params in sorted(output.items()):
        file.write("[%s]\n"%section)
        for line in params:
          file.write(line)
        file.write("\n")
    file.close()
  except:
    print('no prior file available')

