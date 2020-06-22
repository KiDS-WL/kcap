#!/bin/bash 
foldername=/where/the/inifile/is/
ini_file=${foldername}/config/pipeline.ini
multinest_file=/where/the/multinest/file/is/filename # This can in principle be any file with the same columns as the multinest file (doesn't need the weights though). So it can be a list of parameter values with the same ordering that would be in a output multinest file produced by cosmosis usign the same inifile as ${ini_file}. Then code will look for the columns called post and like to find the row with the best fit values and starts from there. Alternatively if this is not given will look in the values.ini file that is named in the ${ini_file} for the start values. 
output_file_name=/where/the/output/will/be/saved/output_name
bestfit_file_name=/where/the/bestfit/values_ini/will/be/saved/bestfit_name_values
bestfit_prior_name=/where/the/bestfit/priors_ini/will/be/saved/bestfit_name_priors # If  --best_fit_value set but this is not set it will still produce somethig with the same name but a "prior_" added to the start

# If -s is removed it will only save the best fit result and not the steps.

python maxlike_cosmosis.py -i ${ini_file} -m ${multinest_file} -o ${output_file_name} -s --best_fit_value ${bestfit_file_name} --best_fit_priors ${bestfit_prior_name} --maxiter 3000  
