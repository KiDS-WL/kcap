#!/bin/sh
# This shouldn't be necessary anymore soon.
export COSMOSIS_SRC_DIR=`python -c "import cosmosis; import os; print(os.path.split(os.path.dirname(cosmosis.__file__))[0])"`