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

file(GLOB SRC_SOLVE counter.f90       mod_blocking_size.f90  omp_bayes_solve.f90     omp_mvp2.f90              omp_sym_solve_ilu.f90   preconditioners.f90  solve_debug.f90
ddl_du.f90        mod_OMP_TOOLS.f90      omp_damax.f90           omp_mvp.f90               omp_sym_solve_ssor.f90  prepare_solve.f90    solve.f90
dense_solve.f90   nag_gen_solve.f90      omp_ddot.f90            OMP_TOOLS.f90             qddot.f90               solve_type.f90
direct_solve.f90  norm_linsys2.f90       omp_gen_solve_diag.f90  par_tools.f90             reduction.f90           ssor_mvp_single.f90
get_dnorm.f90     norm_linsys.f90        omp_gen_solve.f90       omp_preconditioners.f90   p_pos_anz.f90           set_dval.f90         test_matrix.f90
get_norm2.f90     norm_resid.f90         omp_gen_solve_ilu.f90   omp_sym_solve_diag.f90    pre_bicgstab.f90        set_ival.f90         test_symmetry.f90
get_norm.f90      omp_abbruch.f90        omp_gen_solve_ssor.f90  omp_sym_solve.f90         pre_cg.f90              set_lval.f90         test_zero.f90   )
add_library(solve ${SRC_SOLVE})
