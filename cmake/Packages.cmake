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

#packages
if(omp)
   find_package(OpenMP REQUIRED)
   add_compile_options(${OpenMP_Fortran_FLAGS})
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_Fortran_FLAGS}")
   add_definitions(-DfOMP)
endif()
find_package(BLAS REQUIRED)
link_libraries(${BLAS_LIBRARIES})
find_package(LAPACK REQUIRED)
link_libraries(${LAPACK_LIBRARIES})
if(hdf)
   find_package(HDF5 REQUIRED COMPONENTS Fortran Fortran_HL)
   link_libraries(${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})
   include_directories(${HDF5_INCLUDE_DIRS})
   add_definitions(-DHDF)
else()
   add_definitions(-DnoHDF)
endif()
