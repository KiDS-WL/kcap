cmake_minimum_required(VERSION 3.2)

project(kcap VERSION 0.0.1)
    add_subdirectory(cosmosis-standard-library)
    add_subdirectory(HMx)

    add_custom_target(kcap ALL)
    add_dependencies(kcap csl hmx_cosmosis_interface)
    
