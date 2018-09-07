add_custom_target(cl_to_xi_nicaea
                  COMMENT "Running make for cl_to_xi_nicaea."
                  COMMAND ${CMAKE_COMMAND} -E env GSL_INC=${GSL_INCLUDE_DIRS} GSL_LIB=${GSL_LIBRARY_DIRS} FFTW_INCLUDE_DIR=${FFTW_INCLUDES} FFTW_LIBRARY=${FFTW_LIBRARY_DIRS} COSMOSIS_SRC_DIR=${COSMOSIS_SRC_DIR} make
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  DEPENDS GSL FFTW)

add_custom_target(cl_to_xi_nicaea-clean
                  COMMENT "Running make clean for cl_to_xi_nicaea."
                  COMMAND ${CMAKE_COMMAND} -E env COSMOSIS_SRC_DIR=${COSMOSIS_SRC_DIR} make clean
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})