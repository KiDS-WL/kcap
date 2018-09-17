
add_custom_target(camb
                  COMMENT "Running make for CAMB"
                  COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} FC=${FC} MPIFC=${MPIFC} COSMOSIS_SRC_DIR=${COSMOSIS_SRC_DIR} make ${QUIET_MAKE}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                 )
add_custom_target(camb-clean
                 COMMENT "Running make clean for CAMB"
                 COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} FC=${FC} MPIFC=${MPIFC}  COSMOSIS_SRC_DIR=${COSMOSIS_SRC_DIR} make ${QUIET_MAKE} clean
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                )