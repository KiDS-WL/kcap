
add_custom_target(cosmic_emu
                  COMMENT "Running make for cosmic_emu"
                  COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} FC=${FC} MPIFC=${MPIFC} GSL_INC=${GSL_INCLUDE_DIRS} GSL_LIB=${GSL_LIBRARY_DIRS}  COSMOSIS_SRC_DIR=${COSMOSIS_SRC_DIR} make ${QUIET_MAKE}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                 )
add_custom_target(cosmic_emu-clean
                 COMMENT "Running make clean for cosmic_emu"
                 COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} FC=${FC} MPIFC=${MPIFC} GSL_INC=${GSL_INCLUDE_DIRS} GSL_LIB=${GSL_LIBRARY_DIRS}  COSMOSIS_SRC_DIR=${COSMOSIS_SRC_DIR} make ${QUIET_MAKE} clean
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                )