# - Find the CosmoSIS library
# This module finds if CosmoSIS is installed, and sets the following variables
# indicating where it is. Based on FindNumpy.cmake.
#
#  COSMOSIS_FOUND               - was CosmoSIS found
#  COSMOSIS_VERSION             - the version of CosmoSIS found as a string
#  COSMOSIS_VERSION_MAJOR       - the major version number of CosmoSIS
#  COSMOSIS_VERSION_MINOR       - the minor version number of CosmoSIS
#  COSMOSIS_VERSION_PATCH       - the patch version number of CosmoSIS
#  COSMOSIS_VERSION_DECIMAL     - e.g. version 1.6.1 is 10601
#  COSMOSIS_DIR                 - path to CosmoSIS (COSMOSIS_SRC_DIR)
#  COSMOSIS_INCLUDE_DIR         - path to the CosmoSIS include files
#  COSMOSIS_LIB                 - path to the CosmoSIS library

#============================================================================
# Copyright 2012 Continuum Analytics, Inc.
#
# MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
#============================================================================

# Finding NumPy involves calling the Python interpreter
if(CosmoSIS_FIND_REQUIRED)
    find_package(PythonInterp REQUIRED)
else()
    find_package(PythonInterp)
endif()

if(NOT PYTHONINTERP_FOUND)
    set(COSMOSIS_FOUND FALSE)
endif()

execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "import cosmosis as c; import os; print(c.__version__); print(os.path.dirname(c.__file__));"
    RESULT_VARIABLE _COSMOSIS_SEARCH_SUCCESS
    OUTPUT_VARIABLE _COSMOSIS_VALUES
    ERROR_VARIABLE _COSMOSIS_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT _COSMOSIS_SEARCH_SUCCESS MATCHES 0)
    if(CosmoSIS_FIND_REQUIRED)
        message(STATUS "CosmoSIS not found. Try installing with pip install git+https://bitbucket.org/tilmantroester/cosmosis.git@standalone_setup_fixes#egg=cosmosis-standalone")
        message(FATAL_ERROR
                "CosmoSIS import failure:\n${_COSMOSIS_ERROR_VALUE}")
    else()
        message(STATUS "CosmoSIS not found.")   
    endif()
    set(COSMOSIS_FOUND FALSE)
else()
    # Convert the process output into a list
    string(REGEX REPLACE ";" "\\\\;" _COSMOSIS_VALUES ${_COSMOSIS_VALUES})
    string(REGEX REPLACE "\n" ";" _COSMOSIS_VALUES ${_COSMOSIS_VALUES})
    list(GET _COSMOSIS_VALUES 0 COSMOSIS_VERSION)
    list(GET _COSMOSIS_VALUES 1 COSMOSIS_DIR)

    # Make sure all directory separators are '/'
    string(REGEX REPLACE "\\\\" "/" COSMOSIS_DIR ${COSMOSIS_DIR})

    # Get the major and minor version numbers
    string(REGEX REPLACE "\\." ";" _COSMOSIS_VERSION_LIST ${COSMOSIS_VERSION})
    list(GET _COSMOSIS_VERSION_LIST 0 COSMOSIS_VERSION_MAJOR)
    list(GET _COSMOSIS_VERSION_LIST 1 COSMOSIS_VERSION_MINOR)
    list(GET _COSMOSIS_VERSION_LIST 2 COSMOSIS_VERSION_PATCH)
    string(REGEX MATCH "[0-9]*" COSMOSIS_VERSION_PATCH ${COSMOSIS_VERSION_PATCH})
    math(EXPR COSMOSIS_VERSION_DECIMAL
        "(${COSMOSIS_VERSION_MAJOR} * 10000) + (${COSMOSIS_VERSION_MINOR} * 100) + ${COSMOSIS_VERSION_PATCH}")

    find_library(COSMOSIS_LIB
                 NAMES "cosmosis"
                 PATHS "${COSMOSIS_DIR}/datablock/")
    set(COSMOSIS_INCLUDE_DIR "${COSMOSIS_DIR}/datablock/")

    find_package_message(COSMOSIS
        "Found CosmoSIS: version \"${COSMOSIS_VERSION}\" at ${COSMOSIS_DIR}"
        "[${COSMOSIS_DIR}][${COSMOSIS_VERSION}]")
    set(COSMOSIS_FOUND True)
endif()
