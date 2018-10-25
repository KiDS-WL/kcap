include(ExternalProject)

set(FFTWVersion 3.3.7)
set(FFTWMD5 0d5915d7d39b3253c1cc05030d79ac47)

set(FFTW_DIR_HINT ${CMAKE_BINARY_DIR}/extern)
find_package(FFTW 3.0)

# If GSL is not installed, lets go ahead and compile it
if(NOT FFTW_FOUND )
    message(STATUS "FFTW not found, downloading and compiling from source")
    ExternalProject_Add(FFTW
        PREFIX FFTW
        URL http://www.fftw.org/fftw-${FFTWVersion}.tar.gz
        URL_MD5 ${FFTWMD5}
        DOWNLOAD_NO_PROGRESS 1
        CONFIGURE_COMMAND ./configure --prefix=${FFTW_DIR_HINT} 
                            --enable-shared=no --with-pic=yes
                            --enable-silent-rules
                            --quiet
        BUILD_COMMAND           make -j8
        INSTALL_COMMAND         make install
        BUILD_IN_SOURCE 1)
    set(FFTW_USE_STATIC_LIBS TRUE)
    set(FFTW_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/extern/lib/ )
    set(FFTW_INCLUDES ${CMAKE_BINARY_DIR}/extern/include/)
    set(FFTW_LIBRARIES -lfftw3)
else()
    if(NOT TARGET FFTW)
        add_custom_target(FFTW)
    endif()
endif()
