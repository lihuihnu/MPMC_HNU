cmake_minimum_required(VERSION 3.21)
project(MPMC_HNU LANGUAGES CXX)

add_subdirectory(modules/ad)

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests/ad)
endif()
