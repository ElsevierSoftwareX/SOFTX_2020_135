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

if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
   if (NOT debug)
      add_compile_options(-w -O3 -vec_report0 -fpp -axSTPW -fp-model fast=2)
   else()
      add_compile_options(-g -O0 -fpp)
   endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
   if (NOT debug)
      add_compile_options(-fno-second-underscore -march=k8 -O3 -funroll-all-loops -fprefetch-loop-arrays -ffast-math -mno-ieee-fp -DG95 -DCLopt -frepack-arrays -ftree-vectorize -funit-at-a-time -cpp)
   else()
      add_compile_options(-g -O0 -DG95 -cpp)
   endif()
endif()
add_definitions(-DUSE_QDDOT)


#if(NOT use_rm)
#   add_definitions(-DSTBAY)
#endif()
add_definitions(-D${INVTYPE})

if (head_base)
   add_definitions(-Dhead_base)
elseif(pres_base)
   add_definitions(-Dpres_base)
else()
   message(FATAL_ERROR "head_base OR pres_base need to be set")
endif()

add_definitions(-DUSER_${USER} -DPROPS_${PROPS} -DSIMUL_${SIMUL} )

if(NOT plt)
   add_definitions(-DNOPLT)
endif()

if(NOT vtk)
   add_definitions(-DNOVTK)
endif()

