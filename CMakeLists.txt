cmake_minimum_required(VERSION 3.2)

project(HMx LANGUAGES Fortran VERSION 0.0.1)
    set(PROJECT_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
    
    option(BUILD_INTERFACE "Build the CosmoSIS interface" OFF)
    option(BUILD_HMX_EXECUTABLE "Build HMx executable" ON)
    option(BUILD_TESTS "Build tests." ON)

    option(DEBUG "Build with debugging settings" OFF)

    enable_language(Fortran)

    # Adding customized cmake module
    list(APPEND CMAKE_MODULE_PATH  "${PROJECT_SOURCE_DIR}/cmake/Modules/")


    if(BUILD_INTERFACE)
        find_package(CosmoSIS REQUIRED)
    else()
        find_package(CosmoSIS)
    endif()

    set(HMX_FFLAGS  "-Warray-bounds -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none -O3 -fdefault-real-8 -fdefault-double-8")
    if(DEBUG)
        set(HMX_FFLAGS "${HMX_FFLAGS} -Wall -fcheck=all -fbounds-check -fbacktrace -Og")
    endif()
    
    add_library(hmx STATIC  src/constants.f90
                            src/random_numbers.f90
                            src/file_info.f90
                            src/logical_operations.f90
                            src/fix_polynomial.f90
                            src/array_operations.f90
                            src/table_integer.f90
                            src/special_functions.f90
                            src/interpolate.f90
                            src/solve_equations.f90
                            src/string_operations.f90
                            src/calculus_table.f90
                            src/cosmology_functions.f90
                            src/HMx.f90
                            src/Limber.f90
                            src/cosmic_emu_stuff.f90 )

    set(CMAKE_Fortran_FLAGS "${HMX_FFLAGS} -std=gnu -ffree-line-length-none")

    if(COSMOSIS_FOUND)
        add_library(hmx_cosmosis_interface MODULE src/cosmosis_interface.f90)
        add_dependencies(hmx_cosmosis_interface hmx)
        include_directories("${COSMOSIS_INCLUDE_DIR}")
        target_link_libraries(hmx_cosmosis_interface ${COSMOSIS_LIB} hmx)

        set_target_properties(hmx_cosmosis_interface PROPERTIES PREFIX "")
        set_target_properties(hmx_cosmosis_interface PROPERTIES OUTPUT_NAME "cosmosis_interface")

        install(TARGETS hmx_cosmosis_interface
                LIBRARY DESTINATION ${PROJECT_INSTALL_PREFIX}/lib)
    endif()

    if(BUILD_HMX_EXECUTABLE)
        add_executable(hmx_exec src/HMx_driver.f90)
        target_link_libraries(hmx_exec hmx)

        set_target_properties(hmx_exec PROPERTIES OUTPUT_NAME "HMx")

        install(TARGETS hmx_exec
                RUNTIME DESTINATION ${PROJECT_INSTALL_PREFIX}/bin)
    endif(BUILD_HMX_EXECUTABLE)

    if(BUILD_TESTS)
        add_executable(test_gas_gas tests/test_gas_gas.f90)
        add_dependencies(test_gas_gas hmx)
        target_link_libraries(test_gas_gas hmx)

        install(TARGETS test_gas_gas
                RUNTIME DESTINATION ${PROJECT_INSTALL_PREFIX}/bin/tests)
    endif(BUILD_TESTS)

    install(TARGETS hmx
            LIBRARY DESTINATION ${PROJECT_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${PROJECT_INSTALL_PREFIX}/lib)

    # Adds uninstall target
    if(NOT TARGET uninstall)
        configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
                        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                        IMMEDIATE @ONLY)
        add_custom_target(uninstall
                          COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
    endif()
