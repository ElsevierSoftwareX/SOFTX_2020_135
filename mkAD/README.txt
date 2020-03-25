Files in this directory, and their functions:

Makefile		 - creates subdirectories in process_gta{f/p} with independent sets of shemat code for ad of any user/props combination
Makefile.sub		 - to be copied to all of these subdirectories, to call preprocessing(1), ad-tool, postprocessing(2)

tree.csh		 - (1) parse tree of functions, subroutines and modules actually needed, create script to delete all others
preparefiles.csh		 - (1) comment out & delete unhandy & unnecessary constructions (runon-comments, write commands, etc.)
adjust_function_type.py	 - (1) move type declarations from function header to variable declaration (after "implicit none") - prevents tapenade problems
remove_OMP.py		 - (1) transcribe OMP calls into dummy routines
remove_preproc.py	 - (1) transcribe preprocessor statements into dummy routines
remove_dummy.py		 - (2) restore OMP calls and preprocessor statements from dummy routines
use_include_g.py	 - (2) remove unnecessary "g_" and "_ad" from module and include names (for tapenade only)
parser.py		 - (2) find ISIZEs needed for diffsizes.f90 (for tapenade only)

tafdirectives.f90	 - collected taf directives, to be copied to all subdirectories (for taf only)
solve_dummy.f90		 - dummy solve routine to simulate correct in and out variables, but skip ad of solve.f90 (for tapenade only)

f90split.f90.gz		 - to split files into separate functions (Might be needed if TAF -v1 is discontinued, since -v2 creates files with both derivative and re-parsed original code. We want to discard their re-parsed "original" code and continue using ours.)

staf			 - taf script (update as needed)
calc_pres_dummy.f90	 - not needed any more
