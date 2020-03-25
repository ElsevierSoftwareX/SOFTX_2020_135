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

option(omp "OpenMP switch" ON)
option(hdf "HDF5 switch" ON)
option(mpi "MPI switch" OFF)
option(plt "PLT switch" ON) 
option(vtk "VTK switch" ON) 
#option(head "Use head based computation" ON)
#option(pres "Use pressure based computation" OFF)
option(details "Show compilation details" OFF)
option(debug "Compile with debug informations" OFF)
option(use_rm "Use Reverse-Mode in Inversion" OFF)

set(phys_base "head" CACHE INTERNAL "The used physical pressure representation, could be head or pres")
set_property(CACHE phys_base PROPERTY STRINGS "head" "pres")
message(STATUS "phys_base='${phys_base}'")

set(USER "none" CACHE STRING "The Userfunction to be used (standard is none)")
set(PROPS "const" CACHE STRING "The Properties to be used (standard is const)")

set(INVTYPE "STBAY" CACHE STRING "Inversion methods")
set_property(CACHE INVTYPE PROPERTY STRINGS "STBAY" "DSPAC" "PSPAC" "MF_STBAY")
# Solvertype for Nonlinear Coupling
set(NLSOLVETYPE "stdFW" CACHE STRING "Nonlinear Solver")
set_property(CACHE NLSOLVETYPE PROPERTY STRINGS "stdFW" "nwtFW" "nitFW")


set(ADTYPE "g_tap" CACHE STRING "Used AD forward mode directories")
set(AD_RMTYPE "ad_tap" CACHE STRING "Used AD reverse mode directories")

set(SIMUL "sgsim" CACHE STRING "Used Simulation (sgsim)")
if (${phys_base} MATCHES "pres")
   set(pres_base ON)
   set(head_base OFF)
elseif(${phys_base} MATCHES "head")
   set(head_base ON)
   set(pres_base OFF)
else()
   message(FATAL_ERROR "phys_base needs to be set to head OR pres. Add -Dphys_base=head OR -Dphys_base=pres to cmake command [Standard is head]")
endif()

if(use_rm)
   set(matfree ON)
else()
   set(matfree OFF)
endif()

if(details)
   set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

if(debug)
   set(CMAKE_BUILD_TYPE Debug)
endif()
