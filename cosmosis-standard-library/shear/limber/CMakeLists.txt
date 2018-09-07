add_custom_target(limber
                  COMMENT "Running make for limber."
                  COMMAND ${CMAKE_COMMAND} -E env GSL_INC=${GSL_INCLUDE_DIRS} GSL_LIB=${GSL_LIBRARY_DIRS} COSMOSIS_SRC_DIR=${COSMOSIS_SRC_DIR} make
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(limber-clean
                  COMMENT "Running make clean for limber."
                  COMMAND ${CMAKE_COMMAND} -E env GSL_INC=${GSL_INCLUDE_DIRS} GSL_LIB=${GSL_LIBRARY_DIRS} COSMOSIS_SRC_DIR=${COSMOSIS_SRC_DIR} make clean
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})