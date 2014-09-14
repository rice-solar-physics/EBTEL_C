#EBTEL-C

##Authors
EBTEL was originally coded in the Interactive Data Language (IDL) by J.A. Klimchuk, S. Patsourakos, and P.J. Cargill. This version is a translation into the C Programming Language by Will Barnes, Rice University.

For more information regarding the EBTEL model see:

+ <a href="http://adsabs.harvard.edu/abs/2008ApJ...682.1351K">Klimchuk et al. 2008, ApJ, 682:1351-1362</a>
+ <a href="http://adsabs.harvard.edu/abs/2012ApJ...752..161C">Cargill et al. 2012A, ApJ, 752:161</a>
+ <a href="http://adsabs.harvard.edu/abs/2012ApJ...758....5C">Cargill et al. 2012B, ApJ, 758:5</a>

##Model Details
The Enthalpy Based Thermal Evolution of Loops (EBTEL) model allows one to efficiently compute spatially-averaged, time-dependent plasma parameters. It is often desirable to compute solutions for a large number of coronal loops. However, the spatial and temporal scales needed to solve the full _1D_-hydrodynamic equations lead to long computation times. EBTEL allows for quick and accurate solutions to spatially-averaged quantities which allows for an analysis of how the coronal plasma responds to time-dependent heating. Comparisons to _1D_-models have shown (see above publications) that the EBTEL solutions are able to reproduce averaged _1D_-results to a surprising degree of accuracy.

Using a time-dependent heating function, EBTEL computes the plasma response of the temperature, pressure, and number density by solving the _0D_ hydrodynamic equations. These are derived by integrating and energy and mass equations over the coronal and transition region portions of the loop. It is also assumed that the ratios of the average temperature to the apex temperature and the base temperature to the apex temperature are constants.

EBTEL also calculates the differential emission measure as a function of temperature (_DEM(T)_) for the transition region and the corona. Details regarding the specific equations can be found in the above publications.

##Updates in EBTEL-C
The EBTEL-C code offers several advantages over the original IDL code. Perhaps the biggest advantage is the time needed to compute a single run. A typical EBTEL-C run (_~20,000_ s) takes approximately 0.5 seconds as opposed to the _~10_ s that a an EBTEL-IDL run would need (Though it should be noted that this is still efficient compared to a _1D_ hydrodynamic simulation). This drastic increase is speed is due mostly to addition of an adaptive fourth-order Runge-Kutta routine to solve the EBTEL equations. The amount by which the step-size is allowed to vary can be adjusted to allow for greater accuracy or greater speed, whichever is more necessary for the user. 

Additionally, EBTEL-C comes with a more flexible heating function that allows the user to use a variable number of uniform or randomly-occuring heating events or provide an input-file that specifies the heating profile. 

It should be noted that extensive testing has been carried out to ensure that EBTEL-C solutions match those of the EBTEL-IDL solutions such that no additional error is introduced in the translation.

##Downloading and Compiling
The best way to use obtain this code is to clone a copy of this repository on your local machine. If you have `git` installed locally on your machine, to create a working copy type `git clone https://github.com/rice-solar-physics/EBTEL_repo`. Changes may be made periodically to the main EBTEL-C repository. To pull down these changes, but not override any local changes, use `git pull` inside the directory of your working directory.

If the user has a utility installed, the included makefile can be run simply by typing
`make`
at the command line. This creates the object files and links them into the executable 'ebtel' which can be run by typing
`./ebtel`
at the command line. The object files and executable can be cleaned up using the 'clean' option included in the makefile
`make clean`

