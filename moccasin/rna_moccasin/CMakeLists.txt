cmake_minimum_required(VERSION 3.15...3.28)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

# cmake modules for finding external dependencies (none for now)
# set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# add conda environment to search path for finding packages, etc.
if(DEFINED ENV{CONDA_PREFIX})
  set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};$ENV{CONDA_PREFIX})
endif()

# get include/library paths

find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module NumPy REQUIRED)

python_add_library(
    internals
    MODULE
    "include/moccasin/InferParams.cpp"
    "src/internals.cpp"
)
target_include_directories(
    internals
    PUBLIC
    "include/"
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
    "third-party/eigen-3.4.0/"
)
target_compile_features(internals PUBLIC cxx_std_17)
target_compile_definitions(
    internals
    PRIVATE
    NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
    EIGEN_MPL2_ONLY
)

install(
    TARGETS
    internals
    LIBRARY
    DESTINATION
    moccasin/
)
