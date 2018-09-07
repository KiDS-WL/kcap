
add_custom_target(projection
                  COMMENT "Running make for projection"
                  COMMAND ${CMAKE_COMMAND} -E env GSL_INC=${GSL_INCLUDE_DIRS} GSL_LIB=${GSL_LIBRARY_DIRS} COSMOSIS_SRC_DIR=${COSMOSIS_SRC_DIR} make ${QUIET_MAKE}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  DEPENDS GSL)
add_custom_target(projection-clean
                  COMMENT "Running make clean for projection"
                  COMMAND ${CMAKE_COMMAND} -E env GSL_INC=${GSL_INCLUDE_DIRS} GSL_LIB=${GSL_LIBRARY_DIRS} COSMOSIS_SRC_DIR=${COSMOSIS_SRC_DIR} make ${QUIET_MAKE} clean
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
