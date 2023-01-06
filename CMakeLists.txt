
cmake_minimum_required(VERSION 3.20)

project(LC_hypercomplex LANGUAGES CXX)

set(LC_hypercomplex_MAIN_PROJECT OFF)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(LC_hypercomplex_MAIN_PROJECT ON)
endif()

include(FetchContent)

add_library(LC_hypercomplex INTERFACE)

target_include_directories(LC_hypercomplex INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

set_target_properties(LC_hypercomplex PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

FetchContent_Declare(
    LC_math_linalg
    GIT_REPOSITORY https://github.com/LucaCiucci/LC_math_linalg
    GIT_TAG        origin/old_implementation
)
FetchContent_MakeAvailable(LC_math_linalg)

target_link_libraries(LC_hypercomplex INTERFACE LC_math_linalg)


if(LC_hypercomplex_MAIN_PROJECT)
    add_subdirectory(examples)
endif()

