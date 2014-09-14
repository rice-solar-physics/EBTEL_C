#EBTEL-C

##Authors
EBTEL was originally coded in the Interactive Data Language (IDL) by J.A. Klimchuk, S. Patsourakos, and P.J. Cargill. This version is a translation into the C Programming Language by Will Barnes, Rice University.

For more information regarding the EBTEL model see:

+ <a href="http://adsabs.harvard.edu/abs/2008ApJ...682.1351K">Klimchuk et al. 2008, ApJ, 682:1351-1362</a>
+ <a href="http://adsabs.harvard.edu/abs/2012ApJ...752..161C">Cargill et al. 2012A, ApJ, 752:161</a>
+ <a href="http://adsabs.harvard.edu/abs/2012ApJ...758....5C">Cargill et al. 2012B, ApJ, 758:5</a>

##Model Details

##Updates in EBTEL-C

##Compiling
If the user has a utility installed, the included makefile can be run simply by typing
`make`
at the command line. This creates the object files and links them into the executable 'ebtel' which can be run by typing
`./ebtel`
at the command line. The object files and executable can be cleaned up using the 'clean' option included in the makefile
`make clean`

