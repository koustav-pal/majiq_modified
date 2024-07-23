cmake_minimum_required(VERSION 3.15...3.28)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

# cmake modules for finding external dependencies
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# add conda environment to search path for finding packages, etc.
if(DEFINED ENV{CONDA_PREFIX})
  set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};$ENV{CONDA_PREFIX})
endif()


# get include/library paths

find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module NumPy REQUIRED)

set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

find_package(HTSlib REQUIRED)

find_package(ZLIB REQUIRED)
add_subdirectory("third-party/bxzstr-1.2.1")

pybind11_add_module(
    internals
    MODULE
    "include/majiq/internals/py/module.cpp"
    "include/majiq/internals/SJBinsReads.cpp"
    "include/majiq/internals/GFF3.cpp"
)
target_include_directories(
    internals
    PUBLIC
    "include/majiq/"
    "third-party/boost_1_85_0_subset/"
    ${HTSLIB_INCLUDE_DIR}
)
add_dependencies(internals bxzstr)
target_link_libraries(internals PRIVATE ${HTSLIB_LIBRARIES} bxzstr)
target_compile_features(internals PUBLIC cxx_std_17)
target_compile_definitions(
    internals
    PRIVATE
    VERSION_INFO=${SKBUILD_PROJECT_VERSION}
    NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
    )

python_add_library(
    beta_mixture
    MODULE
    "include/majiq/gufuncs/beta_mixture/beta_mixture_module.cpp"
)
target_include_directories(
    beta_mixture
    PUBLIC
    "include/majiq/"
    "third-party/boost_1_85_0_subset/"
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
)
target_compile_features(beta_mixture PUBLIC cxx_std_11)
target_compile_definitions(
    beta_mixture
    PRIVATE
    NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
)

python_add_library(
    _offsets
    MODULE
    "include/majiq/gufuncs/offsets/offsets_module.cpp"
)
target_include_directories(
    _offsets
    PUBLIC
    "include/majiq/"
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
)
target_compile_features(_offsets PUBLIC cxx_std_11)
target_compile_definitions(
    _offsets
    PRIVATE
    NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
)

python_add_library(
    _stats
    MODULE
    "include/majiq/gufuncs/stats/stats_module.cpp"
)
target_include_directories(
    _stats
    PUBLIC
    "include/majiq/"
    "third-party/boost_1_85_0_subset/"
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
)
target_compile_features(_stats PUBLIC cxx_std_11)
target_compile_definitions(
    _stats
    PRIVATE
    NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
)


install(
    TARGETS
    internals
    beta_mixture
    _offsets
    _stats
    LIBRARY
    DESTINATION
    rna_majiq/
)
