#!/bin/bash/
# this prevents CLASS from parallelizing (might be necessary when you run on cluster with queuing system):
#export OMP_NUM_THREADS=1
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
        `# choose the MultiNest sampler (nested sampling)` \
        -m NS \
        `# set an arbitrary but large number of steps (run should converge well before reaching that number!)`
        --NS_max_iter 1000000000 \
        `# flag for using importance nested sampling (more accurate evidence estimates are produced faster)` \
        --NS_importance_nested_sampling True \
        `# for more accurate evidences use 0.3 (for parameter estimation only use 0.8); for consistency with PolyChord this should be set even lower (e.g. 0.01)` \
        --NS_sampling_efficiency 0.3 \
        `# the more live points the smoother the contours, but the runtime scales linearly with this number! A good starting point (inspired by PolyChord settings) is nDims * 25` \
        --NS_n_live_points 350 \
        `# run will finish/is converged if ln(Z_i) - ln(Z_j) <= NS_evidence_tolerance for i>j (0.5 is sufficient for parameter estimation only); for consistency with PolyChord this should be set to 0.01` \
        --NS_evidence_tolerance 0.1
