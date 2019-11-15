cmake_minimum_required(VERSION 3.2)

project(COSEBIS LANGUAGES C CXX VERSION 0.1)
    set(QUIET_MAKE "--no-print-directory")

    # Adding customized cmake module
    list(APPEND CMAKE_MODULE_PATH  "${PROJECT_SOURCE_DIR}/cmake/Modules/")

    find_package(CosmoSIS REQUIRED)

    set(COSMOSIS_SRC_DIR ${COSMOSIS_DIR})

    include(BuildGSL)

    if(APPLE)
        set(LAPACK_LINK "-framework Accelerate")
    else()
        find_package(LAPACK REQUIRED)
	list(GET LAPACK_LIBRARIES 0 LAPACK_LIB)
	get_filename_component(LAPACK_LIB_DIR ${LAPACK_LIB} DIRECTORY)
        set(LAPACK_LINK "-L${LAPACK_LIB_DIR} -llapack -lblas")
	message(STATUS "LAPACK_LINK: ${LAPACK_LIB} ${LAPACK_LIB_DIR}")
    endif()
    message("GSL_INC:" ${GSL_INCLUDE_DIRS})
    message("GSL_LIB:" ${GSL_LIBRARY_DIRS})
    add_custom_target(cosebis
                      COMMENT "Running make for cosebis"
                      COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} GSL_INC=${GSL_INCLUDE_DIRS} GSL_LIB=${GSL_LIBRARY_DIRS} COSMOSIS_SRC_DIR=${COSMOSIS_SRC_DIR} make ${QUIET_MAKE}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                      DEPENDS GSL)
    add_custom_target(cosebis-clean
                      COMMENT "Running make clean for cosebis"
                      COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} GSL_INC=${GSL_INCLUDE_DIRS} GSL_LIB=${GSL_LIBRARY_DIRS} COSMOSIS_SRC_DIR=${COSMOSIS_SRC_DIR} make ${QUIET_MAKE} clean
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
