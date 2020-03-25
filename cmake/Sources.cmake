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

#Module Files - need to be build prior to all other files
file (GLOB SRC_MODULES forward/arrays.f90 forward/mod_*.f90 solve/mod_*.f90 simul/mod_*.f90 hdf5/mod_*.f90) 
file (GLOB SRC_MODULES_AD inverse/mod_*.f90 inverse/${ADTYPE}/g_arrays.f90 ${ADTYPE}/g_mod*.f* props/${PROPS}/${ADTYPE}/g_mod*.f*) 
file (GLOB SRC_MODULES_SM simul/mod_*.f90 simul/enkf/m_*.f*)

#Forward Files
file (GLOB SRC_FORWARD forward/*.f* forward/input/*.f* forward/output/*.f* forward/shemach/*.f* forward/strngut/*.f* forward/ctrlut/*.f*
  forward/temp/*.f* forward/conc/*.f* forward/mathfuncs/*.f*)
file (GLOB SRC_FORWARD_HEAD forward/head/*.f*)
file (GLOB SRC_FORWARD_PRES forward/pres/*.f*)
if (pres_base)
   list(APPEND SRC_FORWARD ${SRC_FORWARD_PRES})
elseif(head_base)
   list(APPEND SRC_FORWARD ${SRC_FORWARD_HEAD})
endif()


#Inverse Files
file(GLOB SRC_INVERSE inverse/*.f* inverse/${ADTYPE}/*.f* ${ADTYPE}/${phys_base}/*.f* blas/${ADTYPE}/*.f*)
if (use_rm)
   file(GLOB SRC_ADHELPER mkAD/ADFirstAidKit/adBuffer.f mkAD/ADFirstAidKit/adStack.c)
   file(GLOB SRC_MODULES_AD_RM inverse/mod_*.f90 inverse/${AD_RMTYPE}/arrays_ad.f90 inverse/${AD_RMTYPE}/mod_*_ad.f* props/${PROPS}/${AD_RMTYPE}/${phys_base}/mod_*_ad.f*)
   list(APPEND SRC_MODULES_AD ${SRC_MODULES_AD_RM})
   file(GLOB SRC_INVERSE_RM inverse/${AD_RMTYPE}/* ${AD_RMTYPE}/${phys_base}/*.f* blas/${AD_RMTYPE}/*.f* props/${PROPS}/${AD_RMTYPE}/${phys_base}/*.f* user/${USER}/${AD_RMTYPE}/${phys_base}/*.f*)
   list(APPEND SRC_INVERSE ${SRC_INVERSE_RM} ${SRC_ADHELPER})
endif()

#HDF FILES
file(GLOB SRC_HDF hdf5/add_cube.f90 hdf5/open_hdf5.f90 hdf5/close_hdf5.f90 hdf5/closeopen_hdf5.f90 hdf5/test_hdf5.f90 hdf5/read_hdf5.f90 hdf5/read_hdf5_int.f90 hdf5/create_hdf5.f90 hdf5/write_all_hdf5.f90 hdf5/add_line.f90   hdf5/write_hdf5.f90 hdf5/write_hdf5_int.f90 hdf5/write_parameter2_hdf5.f90 hdf5/add_plane.f90   hdf5/write_parameter_hdf5.f90  hdf5/read_outt_hdf5.f90 )
#HDF FILES FOR INVERSE
file(GLOB SRC_HDF_AD inverse/mod_inverse.f90 hdf5/write_joutt_hdf5.f90  hdf5/read_joutt_hdf5.f90  hdf5/write_inv_hdf5.f90 )


#Solve Files
file(GLOB SRC_SOLVE solve/counter.f90       solve/mod_blocking_size.f90   solve/omp_mvp2.f90              solve/omp_sym_solve_ilu.f90   solve/preconditioners.f90  solve/solve_debug.f90
solve/ddl_du.f90        solve/mod_OMP_TOOLS.f90      solve/omp_damax.f90           solve/omp_mvp.f90               solve/omp_sym_solve_ssor.f90  solve/prepare_solve.f90    solve/solve.f90
solve/dense_solve.f90   solve/nag_gen_solve.f90      solve/omp_ddot.f90            solve/OMP_TOOLS.f90             solve/qddot.f90               solve/solve_type.f90
solve/direct_solve.f90  solve/norm_linsys2.f90       solve/omp_gen_solve_diag.f90  solve/par_tools.f90             solve/reduction.f90           solve/ssor_mvp_single.f90
solve/get_dnorm.f90     solve/norm_linsys.f90        solve/omp_gen_solve.f90       solve/omp_preconditioners.f90   solve/p_pos_anz.f90           solve/set_dval.f90         solve/test_matrix.f90
solve/get_norm2.f90     solve/norm_resid.f90         solve/omp_gen_solve_ilu.f90   solve/omp_sym_solve_diag.f90    solve/pre_bicgstab.f90        solve/set_ival.f90         solve/test_symmetry.f90
solve/get_norm.f90      solve/omp_abbruch.f90        solve/omp_gen_solve_ssor.f90  solve/omp_sym_solve.f90         solve/pre_cg.f90              solve/set_lval.f90         solve/test_zero.f90 )

file(GLOB SRC_SOLVE_AD   solve/omp_bayes_solve.f90) 

#Simul Files
file(GLOB SRC_SIMUL simul/*.f* simul/enkf/*.f* simul/gs/*.f* simul/${SIMUL}/*.f*)

#User Files
file(GLOB SRC_USER user/${USER}/*.f*)
file(GLOB SRC_USER_AD user/${USER}/${ADTYPE}/${phys_base}/*.f*)

#Props Files
file(GLOB SRC_PROPS props/${PROPS}/*.f*)
list(FILTER SRC_PROPS EXCLUDE REGEX ".*gps\_.*\.f*$")
file(GLOB SRC_PROPS_AD props/${PROPS}/${ADTYPE}/${phys_base}/*.f*)


# Include Directories
include_directories("." "forward" "solve" "simul" "${PROJECT_BINARY_DIR}/generated/" "simul/gs")


#Sources for Forward Build
list(APPEND SRC_FORWARD ${SRC_SOLVE} ${SRC_HDF} ${SRC_PROPS} ${SRC_USER} ${SRC_MODULES})

#Sources for Inverse Build
list(APPEND SRC_INVERSE ${SRC_FORWARD} ${SRC_HDF_AD} ${SRC_MODULES_AD} ${SRC_USER_AD} ${SRC_PROPS_AD} ${SRC_SOLVE_AD})
if(pres_base)
   list(FILTER SRC_INVERSE EXCLUDE REGEX ".*\_unconf.*\.f*$")
endif()
if(NOT ${NLSOLVETYPE} MATCHES "stdFW")
   list(APPEND SRC_FORWARD ${SRC_NONLINEAR} ${SRC_INVERSE})
endif()



list(APPEND SRC_SIMUL ${SRC_FORWARD} ${SRC_MODULES_SM})

file(GLOB read_files */read_*.f* */input/read_*.f* props/*/read_user.f*)
list(FILTER read_files EXCLUDE REGEX ".*forward/input/read_check.f90|.*/input/read_restartFW.f90|.*inverse/read_restartINV.f90")


# External library definitions to be used with e.g. EFCOSS
#BLAS Files
file(GLOB SRC_BLAS blas/*.f*)

#LAPACK Files
file(GLOB SRC_LAPACK lapack/*.f*)


#Directories needed for AD generation
list(APPEND AD_DIRECTORIES forward forward/input forward/output forward/shemach forward/strngut forward/ctrlut
  forward/temp forward/conc forward/mathfuncs hdf5 solve blas)
if(pres)
   list(APPEND AD_DIRECTORIES forward/pres)
else()
   list(APPEND AD_DIRECTORIES forward/head)
endif()

# Define user and property directories for AD generation process (wells3d and basc are required, c.f AD.cmake)
list(APPEND USER_DIRECTORIES wells3d none)# gheexpl)# wells wells6 wells3d_s3w wells3dN_CK wells3dN_CK3B wells3d_fine wells3d_stoch wells3d_stoch_3B wells3d_stoch_s3w gheloop gheexpl)
list(APPEND PROPS_DIRECTORIES basc const bas)# ice kola)# conv basd frac ice kola)

list(APPEND USER_DIRECTORIES_FULL USER_DIRECTORIES)
list(APPEND PROPS_DIRECTORIES_FULL PROPS_DIRECTORIES)
list(TRANSFORM USER_DIRECTORIES_FULL PREPEND "user/")
list(TRANSFORM PROPS_DIRECTORIES_FULL PREPEND "props/")


# All the directories for USER and PROPS variables
file(GLOB ALL_USER_DIRECTORIES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/user user/*)
message(STATUS "Available USER='${ALL_USER_DIRECTORIES}'")
file(GLOB ALL_PROPS_DIRECTORIES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/props props/*)
message(STATUS "Available PROPS='${ALL_PROPS_DIRECTORIES}'")
set_property(CACHE USER PROPERTY STRINGS ${ALL_USER_DIRECTORIES})
set_property(CACHE PROPS PROPERTY STRINGS ${ALL_PROPS_DIRECTORIES})
