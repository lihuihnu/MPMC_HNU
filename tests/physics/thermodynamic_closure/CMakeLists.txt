cmake_minimum_required(VERSION 3.21)
project(mpmc_physics_thermodynamic_closure_tests LANGUAGES CXX)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/physics" physics-library)
enable_testing()

option(MPMC_PHYSICS_CLOSURE_ENABLE_SANITIZERS "Enable ASan/UBSan for supported compilers" OFF)

add_executable(mpmc_physics_closure_tests closure_test.cpp)
target_link_libraries(mpmc_physics_closure_tests PRIVATE mpmc::physics)
set_target_properties(mpmc_physics_closure_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_physics_closure_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_physics_closure_tests PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()
if(MPMC_PHYSICS_CLOSURE_ENABLE_SANITIZERS)
    target_compile_options(mpmc_physics_closure_tests PRIVATE -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_physics_closure_tests PRIVATE -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

foreach(case binary_reference ternary_reference unavailable_policy indeterminate_policy)
    add_test(NAME physics.thermodynamic_closure.${case}
             COMMAND mpmc_physics_closure_tests ${case})
    set_tests_properties(physics.thermodynamic_closure.${case} PROPERTIES
        LABELS "physics;thermodynamic_closure" TIMEOUT 180)
endforeach()
