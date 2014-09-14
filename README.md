#EBTEL-C

##AUTHORS
EBTEL was originally coded in the Interactive Data Language (IDL) by J.A. Klimchuk, S. Patsourakos, and P.J. Cargill. This version is a translation into the C Programming Language.

For more information regarding the EBTEL model see:

+<a href="http://adsabs.harvard.edu/abs/2008ApJ...682.1351K">Klimchuk et al. 2008, ApJ, 682:1351-1362</a>

+<a href="http://adsabs.harvard.edu/abs/2012ApJ...752..161C">Cargill et al. 2012A, ApJ, 752:161</a>

+<a href="http://adsabs.harvard.edu/abs/2012ApJ...758....5C">Cargill et al. 2012B, ApJ, 758:5</a>

COMPILING
------------------
If the user has a 'make' utility installed, the included makefile can be run simply by typing
$ make
at the command line. This creates the object files and links them into the executable 'ebtel' which can be run by typing
$ ./ebtel
at the command line. The object files and executable can be cleaned up using the 'clean' option included in the makefile
$ make clean

DOWNLOADING
------------------

DOCUMENTATION
------------------
This directory contains 6 pieces of source code and a text file which drives the executable. Below is a list of each file and their respective functions.

-ebtel_functions_loop.c
--ebtel_loop_solver
--ebtel_kpar_set
--ebtel_rad_loss
--ebtel_calc_tr_dem
--ebtel_heating

-ebtel_functions_param.c
--ebtel_calc_c1
--ebtel_calc_c2
--ebtel_calc_c3
--ebtel_calc_lambda
--ebtel_calc_abundance
--ebtel_calc_ic

-ebtel_functions_solvers.c
--ebtel_euler
--ebtel_rk
--ebtel_rk_adapt
--ebtel_rk_derivs

-ebtel_functions_util.c
--ebtel_print_header
--ebtel_file_writer
--ebtel_linspace
--ebtel_colon_operator
--ebtel_avg_val
--ebtel_max_val
--ebtel_min_val
--ebtel_free_mem

-ebtel_functions.h

-ebtel_main.c

The code is reasonably well documented. See the above publications for more information.
