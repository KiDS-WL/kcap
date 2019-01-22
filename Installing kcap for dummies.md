Here is a guide to what you need to do to install kcap and what issues you might encounter. 
Setting a `conda` environment is explained at the end of the document.

Choose the location you want the kcap folder to sit and clone it:

```
git clone https://github.com/KiDS-WL/kcap.git
```
    
Enter the kcap folder and follow the installation guide in the readme file and use `python build.py` command to build and install all the relevant modules. But that might not work! If it doesn’t these are things to try:

Set the paths to CC, CXX (and MPIFC if you need parallel run support using mpi), try this in bash shell:

```
export CC=gcc

export CXX= g++

export MPIFC=mpif90
```

In mac specially you might not have the latest gcc compiler. You can use homebrew to get the latest:

```
brew install gcc
```

Take a note of where gcc is installed. In mac you might need to specify the version of gcc and the path to it, for example for me the paths look like this:
```
export CC=/usr/local/Cellar/gcc/8.2.0/bin/gcc-8

export CXX= /usr/local/Cellar/gcc/8.2.0/bin/g++-8
```

After setting these try `python build.py` again. If this doesn’t work first uninstall it with:

```
python build.py --clean
```
And try again. Then test kcap using:
```
cosmosis example/example_b.ini
```
If you get an error about cosmosis, then first uninstall cosmosis:

```
pip uninstall cosmosis_standalone
```

And then build kcap again

```
python build.py
```

Run cosmosis again, hopefully it works!

Now for mpi you’ll probably need to install mpi4py:

```
pip install mpi4py
```

Test mpi with: 

```
mpirun -n 4 cosmosis --mpi examples/example_b.ini
```

## Setting a conda environment:

Get miniconda from here: https://conda.io/miniconda.html (This link is broken!)
Installation guide also here: https://conda.io/docs/user-guide/install/linux.html (This link is broken!)

Download the miniconda to the machine that you want to install it on. Note if you download this on a different machine and then copy it to the destination you might get an error similar to this: 

ERROR: size of Miniconda3-latest-Linux-x86_64.sh should be     69826864 bytes

So instead find the link to the miniconda version that is applicable to your machine and use this command:
```
wget -P /address/to/the/Download/Folder/ "full link to the file to be downloaded"
```

Then following the instructions on  https://conda.io/docs/user-guide/install/linux.html type:
```
bash ./Name_of_the_miniconda_file.sh
```

Create an environment (my_env) and set python version to 3:
```
conda create -n my_env python=3 
```

Activate the environment:
```
conda activate my_env
```

You can activate this environment later on and use it with kcap. 

To deactivate write:

```
conda deactivate
```
