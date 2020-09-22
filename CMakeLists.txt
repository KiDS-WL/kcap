cmake_minimum_required(VERSION 3.2)
cmake_policy(SET CMP0002 OLD)

project(kcap VERSION 0.0.1)
    add_subdirectory(cosmosis-standard-library)
    add_subdirectory(cosebis)

    add_custom_target(kcap ALL)
    add_dependencies(kcap csl cosebis)
    
