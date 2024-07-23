# FindHTSlib.cmake
#
# Adapted from FindBlosc.cmake <https://github.com/Intel-HLS/TileDB/blob/master/cmake/Modules/FindBLOSC.cmake>
#
# The MIT License
#
# Copyright (c) 2016 MIT and Intel Corporation
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Finds the HTSLib library. This module defines:
#   - HTSLIB_INCLUDE_DIR, directory containing headers
#   - HTSLIB_LIBRARIES, the Blosc library path
#   - HTSLIB_FOUND, whether Blosc has been found

# Find header files
if(DEFINED ENV{HTSLIB_INCLUDE_DIR})
  find_path(
      HTSLIB_INCLUDE_DIR htslib/hts.h
      PATHS $ENV{HTSLIB_INCLUDE_DIR}
      NO_DEFAULT_PATH
  )
elseif(HTSLIB_SEARCH_HEADER_PATHS)
  find_path(
      HTSLIB_INCLUDE_DIR htslib/hts.h
      PATHS ${HTSLIB_SEARCH_HEADER_PATHS}
      NO_DEFAULT_PATH
  )
else()
  find_path(HTSLIB_INCLUDE_DIR htslib/hts.h)
endif()

# Find library
if(DEFINED ENV{HTSLIB_LIBRARY_DIR})
  find_library(
      HTSLIB_LIBRARIES NAMES hts
      PATHS $ENV{HTSLIB_LIBRARY_DIR}
      NO_DEFAULT_PATH
  )
elseif(HTSLIB_SEARCH_LIB_PATH)
  find_library(
      HTSLIB_LIBRARIES NAMES hts
      PATHS ${HTSLIB_SEARCH_LIB_PATH}
      NO_DEFAULT_PATH
  )
else()
  find_library(HTSLIB_LIBRARIES NAMES hts)
endif()

if(HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARIES)
  message(STATUS "Found HTSLib: ${HTSLIB_LIBRARIES}")
  set(HTSLIB_FOUND TRUE)
else()
  set(HTSLIB_FOUND FALSE)
endif()

if(HTSLIB_FIND_REQUIRED AND NOT HTSLIB_FOUND)
  message(FATAL_ERROR "Could not find the HTSLib library.")
endif()

#include(FindPackageHandleStandardArgs)
#find_package_handle_standard_args(HTSLIB DEFAULT_MSG)
