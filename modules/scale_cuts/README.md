README for load_cterm & scale_cuts modules for KCAP
===================================================
Chieh-An Lin  
2019.11.15


Description
-----------

Ta da!

This contains 
1. a module to read in c-terms and to stock them in the data block (which will be used in the theory pipeline) and 
2. a module to do scale cuts, in which we also define which statistics to include in the likelihood.


Requirement
-----------

The scale cuts module might not work without the newest version of the theory pipeline. 
An up-to-date version can be found on the `bandpower_likelihood` branch of Marika's 
repository: https://bitbucket.org/marika_a/cosebis_cosmosis/src/bandpower_likelihood/


Default test
------------

The `ptest.ini` should provide a default test using KV450 data. The scale cuts are 
set to KV450-like: keeping 7 smallest scales in xi_+ & 6 largest scales in xi_-.

Cuts are defined in the `[scale_cuts]` section. More details are set up in a separate
file `scale_cuts.ini`. Check the instruction in that file.

