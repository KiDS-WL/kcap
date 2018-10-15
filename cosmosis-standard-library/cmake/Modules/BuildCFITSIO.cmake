include(ExternalProject)

set(CFITSIOVersion 3.450)
set(CFITSIOSHA256 bf6012dbe668ecb22c399c4b7b2814557ee282c74a7d5dc704eb17c30d9fb92e)

string(REPLACE "." "" CFITSIOVersionFilename ${CFITSIOVersion})

set(CFITSIO_DIR_HINT ${CMAKE_BINARY_DIR}/extern)
find_package(CFITSIO 3.2)

# If GSL is not installed, lets go ahead and compile it
if(NOT CFITSIO_FOUND )
    message(STATUS "CFITSIO not found, downloading and compiling from source")
    ExternalProject_Add(CFITSIO
                        PREFIX CFITSIO
                        URL https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio${CFITSIOVersionFilename}.tar.gz
                        URL_HASH  SHA256=${CFITSIOSHA256}
                        DOWNLOAD_NO_PROGRESS 1
                        CONFIGURE_COMMAND   ./configure
                                            --prefix=${CFITSIO_DIR_HINT}
                                            --enable-silent-rules
                                            --quiet
                        BUILD_COMMAND           make -j8 > out.log 2>&1
                        INSTALL_COMMAND         make install
                        BUILD_IN_SOURCE 1
                       )

    set(CFITSIO_LIBRARY_DIR ${CMAKE_BINARY_DIR}/extern/lib/ )
    set(CFITSIO_INCLUDE_DIR ${CMAKE_BINARY_DIR}/extern/include/)
    set(CFITSIO_LIB -lcfitsio)
    set(CFITSIO_FOUND True)
else()
    if(NOT TARGET CFITSIO)
        add_custom_target(CFITSIO)
    endif()
endif()
