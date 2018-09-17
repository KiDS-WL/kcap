cmake_minimum_required(VERSION 3.2)

project(CSL LANGUAGES C CXX Fortran VERSION 0.0.8)
    set(FC "gfortran" CACHE STRING "Fortran compiler.")

    set(QUIET_MAKE "--no-print-directory")

    add_custom_target(csl ALL)

    # Adding customized cmake module
    list(APPEND CMAKE_MODULE_PATH  "${PROJECT_SOURCE_DIR}/cmake/Modules/")

    find_package(CosmoSIS REQUIRED)

    set(COSMOSIS_SRC_DIR ${COSMOSIS_DIR})

    include(BuildGSL)
    include(BuildFFTW)
    include(BuildCFITSIO)

    add_subdirectory(boltzmann/camb)
    add_subdirectory(boltzmann/halofit)
    add_subdirectory(boltzmann/halofit_takahashi)

    add_subdirectory(shear/cl_to_xi_nicaea)

    add_subdirectory(structure/projection)

    add_dependencies(csl camb halofit cl_to_xi_nicaea projection)

    add_custom_target(csl-clean)
    add_dependencies(csl-clean camb-clean halofit-clean cl_to_xi_nicaea-clean projection-clean)

    
