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
The EBTEL-C code offers several advantages over the original IDL code. Perhaps the biggest advantage is the time needed to compute a single run. A typical EBTEL-C run (_~20,000_ s) takes approximately 0.5 seconds as opposed to the _~10_ s that a an EBTEL-IDL run would need (Though it should be noted that this is still efficient compared to a _1D_ hydrodynamic simulation). This drastic increase in speed is due mostly to addition of an adaptive fourth-order Runge-Kutta routine to solve the EBTEL equations. The amount by which the step-size is allowed to vary can be adjusted to allow for greater accuracy or greater speed, whichever is more necessary for the user. 

Additionally, EBTEL-C comes with a more flexible heating function that allows the user to use a variable number of uniform or randomly-occuring heating events or provide an input-file that specifies the heating profile. 

It should be noted that extensive testing has been carried out to ensure that EBTEL-C solutions match those of the EBTEL-IDL solutions such that no additional error is introduced in the translation.

![Example EBTEL-C run showing resulting temperature and density profiles from an impulsive heating event](ebtelC_example.png)

##Dependencies
EBTEL-C uses an XML configuration file system to minimize errors that result from poorly formatted input files and allow for easier input readability. EBTEL-C uses the XML C parser toolkit `libxml2` which provides a number of useful functions and datatypes for parsing structured XML files. The toolkit, which is essentially a collection of header files, can be obtained <a href="http://xmlsoft.org/downloads.html">here.</a> Additionally, Mac users can obtain the library using the MacPorts package manager by using `sudo port install libxml2`. Linux users can obtain the library via the built-in Aptitude package manager.

**NOTE: If you use the build procedure described below to compile EBTEL-C, you must first change the include location of the `libxml2` library in `src/makefile`.** At the top of the file, the variable `IFLAGS` specifies the location of your local copy of the `libxml2` directory containing all of the necessary header files. For example, if the `libxml2` directory is in `/opt/local/include`, then the line in your makefile would be:

+ `IFLAGS=-I /opt/local/include/libxml2`

You of course may also choose to compile EBTEL-C by hand or write your own `makefile`.

##Downloading and Compiling
Linux and Mac users should be able to compile and run EBTEL-C in the terminal. Windows users should compile and run the code in the Cygwin environment (<a href="https://www.cygwin.com/">available here</a>). The best way to obtain this code is to clone a copy of this repository on your local machine. If you have `git` installed locally, create a working copy by typing `git clone https://github.com/rice-solar-physics/EBTEL_C.git` at the command line. Changes may be made periodically to the main EBTEL-C repository. To pull down these changes, but not override any local changes, use `git pull` inside of your working directory. You can also simply download a compressed file containing all of the source code if you do not wish to bother with the version control.

To compile EBTEL-C, switch to the `build` directory and run `./build`. This uses a makefile in `src` to compile the source code and place an executable called `ebtel` in the `bin` directory. Additionally, running './clean' in `build` removes the executable and all of the object files created at compile time.

To run EBTEL-C, simply run `./ebtel` <b>in the `bin` directory.</b> EBTEL-C also accepts two optional command line arguments: (1) `quiet` which silences the header printed by default and (2) a custom configuration filename. If no filename is specified, the default filename `../config/ebtel_config.xml` will be used. Custom filename paths should all be relative to the EBTEL-C root directory. The order of the two arguments is arbitrary. All of the following are valid calls of the EBTEL-C executable:

+ `./ebtel`
+ `./ebtel quiet`
+ `./ebtel quiet ../config/new_ebtel_config_file.xml`
+ `./ebtel ../config/new_ebtel_config_file.xml`  

##Configuring Input Parameters
As stated above, EBTEL-C uses an XML configuration file system as opposed to a traditional text file input configuration. XML files allow for increased readability because the context of each parameter (the so-called "node name") is included in the input configuration. Additionally, the order of the input parameters in the configuration file is arbitrary since XML uses keywords rather than the order of the values to associate specific values with specific tags. A sample configuration file is provided in `config/ebtel_config.xml`. Below is a list of the input parameters set by the input configuration file along with a brief description. 

+ General input parameters
  + total_time (s) -- the total amount of time allotted for the simulation
  + tau (s) -- time step; static for Euler and Runge-Kutta solvers; starting time step for adaptive Runge-Kutta solver
  + loop_length (Mm) -- loop half-length; measured from the base of the transition region to the loop apex
  + usage_option -- `dem` include DEM calculation, `no_dem` leave out DEM calculation, `nt_ebeam` include non-thermal electron heating term in pressure equation (this option has NOT been extensively tested with EBTEL-C), `rad_ratio` compute radiation and heat flux ratios and perform DEM calculation. Note that the first and fourth options require longer compute times. If you are only interested in the temperature and density profiles, the second option (`no_dem`) is the recommended choice.
  + rad_option -- choose how to calculate the radiative loss function: `rk` use Raymond-Klimchuk loss function or `rtv` use Rosner-Tucker-Vaiana loss function
  + dem_option -- method for calculating transition region DEM; use either the `new` or `old` method; see <a href="http://adsabs.harvard.edu/abs/2008ApJ...682.1351K">Klimchuk et al. (2008)</a> for details on the differences between these two options
  + heat_flux_option -- EBTEL uses the Spitzer-Harm formula to calculate the heat flux; a flux limiter can be applied by using the `dynamic` option; to calculate the heat flux without the flux limit, use the `classical` option.
  + solver -- three different solvers are available in EBTEL-C: `euler` Euler solver, `rk4` 4th-order Runge-Kutta solver, `rka4` adaptive method coupled to 4th-order Runge-Kutta solver.
  + ic_mode -- choose how to calculate the initial conditions: `st_eq` static equilibrium calculation, `force` initial conditions read in from input file, `scaling` scaling laws calculation.
  + rka_error -- value that defines the allowed error tolerance in the adaptive time step routine; 1.0e-6 is the recommended value
  + index_dem -- index that defines temperature range over which the transition region DEM is computed; 451 is the recommended value.
  + T0 (K) -- temperature at time _t=0_ s
  + n0 (cm^-3) -- number density at time _t=0_ s
+ Heating parameters
  + heating_shape -- shape of the heating pulse; three possible options: `triangle`, `square`, or `gaussian`
  + num_events -- number of heating events (Note that The actual number of events is given by the number of start times that fall within [0,total time])
  + t_start (s) -- time at which the first heating event begins
  + t_pulse_half (s) -- duration of heating event divided by two; for a triangular pulse, this is the time between the start of the event and when the heating amplitude reaches its maximum; for a square pulse, this is just half of the event duration; for a gaussian pulse. this is the sigma parameter.
  + h_back (erg cm^-3 s^-1) -- background heating value
  + t_start_switch -- `uniform`: gives `num_events` heating events beginning at `t_start` separated by 2*`t_pulse_half`; `normal`: selects N start times from a normal distribution with the given mean `mean_t_start` and standard deviation `std_t_start`; `file`: read in start times from values in the `start_time_array` node
    + mean_t_start (s) -- mean value of normally distributed heating event start times
    + std_t_start (s) -- standard deviation of normally distributed start times
  + amp_switch -- `uniform`: gives `num` heating events with uniform amplitude `h_nano`; `power_law`: selects `num_events` heating amplitudes from a power law distribution with index `alpha` and bounds [`amp0`,`amp1`]; `file`: read in amplitudes from values in the `amp_array` node.
    + h_nano (erg cm^-3 s^-1) -- maximum heating amplitude.
    + alpha -- power law index for heating event amplitude distribution
    + amp0 (erg cm^-3 s^-1) -- lower bound on amplitude power law distribution
    + amp1 (erg cm^-3 s^-1) -- upper bound on amplitude power law distribution
  + t_end_switch -- `uniform`: computes end time by adding 2*`t_pulse_half` to each start time giving uniform width to each event; `file`: reads in end times from values in the `end_time_array` node.

##Reporting Bugs and Issues
If you find any bugs or have any concerns about the code, create an Issue or submit a pull request. Questions can also be directed to `will (dot) t (dot) barnes (at) rice (dot) edu`.
