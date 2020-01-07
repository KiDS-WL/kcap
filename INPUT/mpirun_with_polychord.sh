#!/bin/bash/
# this prevents CLASS from parallelizing (recommended for PolyChord):
export OMP_NUM_THREADS=1
# run with MPI:
# if you don't want to / can't run with MPI support, just remove everything before 'python'!
# set number of cores requested for MPI run:
ncores=32
# set absolute path to the root folder of Monte Python (usually 'montepython_public':
mpirun -np $ncores python /your/path/to/montepython_public/montepython/MontePython.py run \
        `# supply relative path from working directory (wd) to param-file` \
	-p your/path/to/K1K_COSEBIs.param \
	`# supply relative path from wd to output folder` \
        -o your/path/to/the/output/folder/ \
	`# supply relative path from wd to correctly set config-file (otherwise default.conf from MontePython will be used)` \
        --conf /your/path/to/your_config.conf \
	`# choose the PolyChord sampler (nested sampling) with its default settings, i.e. --PC_nlive = nDims*25 and --PC_precision_criterion = 0.001` \
        -m PC
