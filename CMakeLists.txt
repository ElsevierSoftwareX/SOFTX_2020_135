# MIT License
#
# Copyright (c) 2020 SHEMAT-Suite
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(SHEMAT-Suite VERSION 9.00 DESCRIPTION "Freely" LANGUAGES Fortran)

enable_language(Fortran C CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake/)

# Executable sources
file(GLOB SRC_MAIN_FW shem_fw.f90)
file(GLOB SRC_MAIN_AD shem_ad.f90)
file(GLOB SRC_MAIN_DD shem_dd.f90)
file(GLOB SRC_MAIN_SM shem_sm.f90)


include(cmake/Options.cmake)
include(cmake/Sources.cmake)
include(cmake/Flags.cmake)
include(cmake/Packages.cmake)
include(cmake/Dependencies.cmake)
include(cmake/AD.cmake)



set(FW_NAME "shem_fw_${PROPS}.x")
add_executable(${FW_NAME} ${SRC_MAIN_FW} ${SRC_FORWARD})

set(AD_NAME "shem_ad_${PROPS}.x")
add_executable(${AD_NAME} EXCLUDE_FROM_ALL ${SRC_MAIN_AD} ${SRC_INVERSE})

set(DD_NAME "shem_dd_${PROPS}.x")
add_executable(${DD_NAME} EXCLUDE_FROM_ALL ${SRC_MAIN_DD} ${SRC_INVERSE})

set(SM_NAME "shem_sm_${PROPS}.x")
add_executable(${SM_NAME} EXCLUDE_FROM_ALL ${SRC_MAIN_SM} ${SRC_SIMUL})

# Library Sources
file(GLOB SRC_LIB_SHEMAT shem_fl.f90)
file(GLOB SRC_LIB_SHEMAT_PARAMETER shem_parameter.f90)

#Shortcut Make Targets
add_custom_target(fw DEPENDS ${FW_NAME})
add_custom_target(ad DEPENDS ${AD_NAME})
add_custom_target(dd DEPENDS ${DD_NAME})
add_custom_target(sm DEPENDS ${SM_NAME})

target_compile_definitions(${AD_NAME} PUBLIC -DAD -Dset${ADTYPE})
if (NOT ${NLSOLVETYPE} MATCHES "stdFW") 
   target_compile_definitions(${FW_NAME} PUBLIC -D${NLSOLVETYPE} -DAD -Dset${ADTYPE})
   target_include_directories(${FW_NAME} PUBLIC nonlinear/nitsol/Nitsol)
else()
   add_definitions(-DstdFW)
endif()


if(use_rm)
   target_compile_definitions(${AD_NAME} PUBLIC -DJACOBI_FREE -DAD_RM -Dset${AD_RMTYPE})
endif()
target_compile_definitions(${DD_NAME} PUBLIC -DAD -Dset${ADTYPE})


add_library("shemat_oed" EXCLUDE_FROM_ALL ${SRC_LIB_SHEMAT} ${SRC_INVERSE} ${SRC_BLAS} ${SRC_LAPACK})
target_compile_definitions("shemat_oed" PUBLIC -DAD -Dset${ADTYPE})
target_compile_options("shemat_oed" PUBLIC -fPIC)


add_library("shemat_parameter" EXCLUDE_FROM_ALL ${SRC_LIB_SHEMAT_PARAMETER} ${SRC_INVERSE} ${SRC_BLAS} ${SRC_LAPACK})
target_compile_definitions("shemat_parameter" PUBLIC -DMF_STBAY -DAD -Dset${ADTYPE} -Dmatvec -DAD_RM -DJACOBI_FREE)
target_compile_options("shemat_parameter" PUBLIC -fPIC)


add_custom_target(lib_oed DEPENDS "shemat_oed")
add_custom_target(lib_parameter DEPENDS "shemat_parameter")
