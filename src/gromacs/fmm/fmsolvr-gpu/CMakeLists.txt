#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

if (NOT DEFINED FMM_STANDALONE)

if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # not using regular Clang or AppleClang
    ADD_DEFINITIONS ( -fabi-version=0 )
endif()

file(GLOB GPUFMM_SOURCES *.cu)
list(REMOVE_ITEM GPUFMM_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/fmmtest.cu")
set(GPUFMM_SOURCES ${GPUFMM_SOURCES} PARENT_SCOPE)

else()

cmake_minimum_required(VERSION 3.2)

# project name
project (fmsolvrv2)
find_package(CUDA "10.0")
SET(CUDA_SEPARABLE_COMPILATION "ON" )
#find_cuda_helper_libs(cudadevrt)

add_definitions(-DCUDADEBUG)
#add_definitions(-DGMX_FMM_DOUBLE)

# setup include directories
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}")
CUDA_INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}")
SET(CMAKE_CXX_FLAGS)
string(APPEND CMAKE_CXX_FLAGS "--std=c++11 -msse3 -mavx -O3 -fopenmp -ffast-math -funroll-loops -fabi-version=6")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3")
#list(APPEND CMAKE_CXX_FLAGS "-O3")
#SET(CUDA_PROPAGATE_HOST_FLAGS "OFF")
#SET(CMAKE_BUILD_TYPE "Debug")

set(NVCC_FLAGS)
#list(APPEND NVCC_FLAGS "-lineinfo")
#list(APPEND NVCC_FLAGS "--maxrregcount=64")
list(APPEND NVCC_FLAGS "-O3")
list(APPEND NVCC_FLAGS "-Wno-deprecated-declarations")
#list(APPEND NVCC_FLAGS -Xcompiler -fPIC)
#list(APPEND NVCC_FLAGS -Xcompiler -fopenmp)
#list(APPEND NVCC_FLAGS "-gencode;arch=compute_35,code=sm_35")
#list(APPEND NVCC_FLAGS "-gencode;arch=compute_52,code=sm_52")
list(APPEND NVCC_FLAGS "-gencode;arch=compute_61,code=sm_61")
list(APPEND NVCC_FLAGS "-gencode;arch=compute_70,code=sm_70")
list(APPEND NVCC_FLAGS "-gencode;arch=compute_75,code=sm_75")
list(APPEND NVCC_FLAGS "-gencode;arch=compute_75,code=compute_75")
#list(APPEND NVCC_FLAGS "-g;-G")
SET(CUDA_NVCC_FLAGS "${NVCC_FLAGS}")

file(GLOB FMM_SOURCES *.cu)
list(REMOVE_ITEM FMM_SOURCES "${PROJECT_SOURCE_DIR}/fmmtest.cu")
list(REMOVE_ITEM FMM_SOURCES "${PROJECT_SOURCE_DIR}/fmsolvr.cu")

#LINK_DIRECTORIES(/home/bkohnke/fmsolvr-build-gromacsv2 /usr/lib64)
ADD_LIBRARY(testdata testdata.cpp)
#CUDA_ADD_LIBRARY(libfmm SHARED ${FMM_SOURCES})
CUDA_ADD_LIBRARY(libfmm ${FMM_SOURCES})
CUDA_ADD_EXECUTABLE(cuda_fmsolvr fmmtest.cu)
#TARGET_LINK_LIBRARIES(cuda_fmsolvr testdata libfmm gomp ${CUDA_cudadevrt_LIBRARY})
TARGET_LINK_LIBRARIES(cuda_fmsolvr testdata libfmm gomp)





endif()
