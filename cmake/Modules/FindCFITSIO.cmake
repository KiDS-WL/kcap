# - Find the CFITSIO library
#
# Usage:
#   find_package(CFITSIO [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   CFITSIO_FOUND               ... true if cfitsio is found on the system
#   CFITSIO_LIB                 ... full path to cfitsio library
#   CFITSIO_INCLUDE_DIR         ... cfitsio include directory


include(FindPkgConfig)
# Check if we can use PkgConfig
find_package(PkgConfig)
#Determine from PKG
if(PKG_CONFIG_FOUND AND NOT CFITSIO_ROOT)
    pkg_check_modules( PKG_CFITSIO QUIET "cfitsio")
endif()

if(CFITSIO_ROOT_DIR)
    #find libs
    find_library(CFITSIO_LIB
                 NAMES "cfitsio"
                 PATHS ${CFITSIO_ROOT_DIR}
                 PATH_SUFFIXES "lib"
                 NO_DEFAULT_PATH
                )
    #find includes
    find_path(CFITSIO_INCLUDE_DIR
              NAMES "fitsio.h"
              PATHS ${CFITSIO_ROOT_DIR}
              PATH_SUFFIXES "include"
              NO_DEFAULT_PATH
             )
else()
    find_library(CFITSIO_LIB
                 NAMES "cfitsio"
                 PATHS ${PKG_CFITSIO_LIBRARY_DIRS} ${LIB_INSTALL_DIR} ${CFITSIO_DIR_HINT}
                )

    find_path(CFITSIO_INCLUDE_DIR
              NAMES "fitsio.h"
              PATHS ${PKG_CFITSIO_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR} ${CFITSIO_DIR_HINT}
             )
endif(CFITSIO_ROOT_DIR)

get_filename_component(CFITSIO_LIBRARY_DIR ${CFITSIO_LIB} DIRECTORY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CFITSIO DEFAULT_MSG
                                  CFITSIO_INCLUDE_DIR CFITSIO_LIB)
mark_as_advanced(CFITSIO_INCLUDE_DIR CFITSIO_LIB)
