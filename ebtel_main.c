/************************************************************************************

FILENAME: ebtel_main.c

AUTHOR Will Barnes

DATE: created: 7 March 2014

DESCRIPTION: This file makes the function call to the primary function in ebtel_functions.c.
It first defines some basic parameters of the EBTEL model including the input parameters as well
as the constants used throughout. Following the necessary function call, it saves the data.

-------------------------------------------------------------------------------------
EBTEL (Enthalpy Based Thermal Evolution of Loops) computes 0-D
hydrodynamic equations.  This software was originally developed by Klimchuk et al,. (2008) and
later improved upon by Cargill et al., (2012a,b). It solves the set of hydrostatic equations by 
averaging over the loop half-length and then integrating in time. Additional details are available
in the two references given above.

Note on variable correspondence with Klimchuk et al. (2008)
   r1 = c_3
   r2 = c_2
   r3 = c_1
   f, ff = F_0
   f_eq, ff_eq = - R_tr
   dem_eq = DEM_se

INTENSITIES:
   For observations in which temperature response function, G(T), has units of 
   DN s^-1 pix^-1 cm^5 and the loop diameter, d, is larger than the pixel dimension, l_pix:

      I_cor_perp = d/(2L)* Int{G(T)*dem_cor(T)*dT}
      I_tr_perp = d/l_pix * Int{G(T)*dem_tr(T)*dT}
      I_tr_parallel = Int{G(T)*dem_tr(T)*dT} ,

   for lines-of-sight perpendicular and parallel to the loop axis.  I_tr_perp assumes that 
   the transition region is thinner than l_pix.

MISCELLANEOUS COMMENTS:
   Runs much more quickly if the transition region DEM is not computed.
   Speed can be increased by increasing the minimum DEM temperature from 10^4 to, say, 10^5 K 
      or by decreasing the maximum DEM temperature from 10^8.5 to, say, 10^7.5 
      (search on 450 and 451).
   The equilibrium base heat flux coefficient of 2/7 is appropriate for uniform heating;
      a coefficient of 4/7 is more appropriate for apex heating.
   To have equal amounts of thermal and nonthermal heating:  flux_nt = heat*length.
   It is desirable to have a low-level background heating during the cooling phase so that the 
      coronal temperature does not drop below values at which the corona DEM is invalid.
   r1 = c_3 = 0.7 gives more accurate coronal evolution than the original 0.5, especially in 
      the late phase of cooling.  However, it produces excess DEM at the very hottest temperatures 
      during impulsive events, and the transition region DEM is somewhat elevated.  We have 
      therefore introduced r1_tr = 0.5, which provides a more accurate transition region DEM at 
      the same time that r1 = 0.7 provides a more accurate radiative cooling.
   v = (c_3/c_2)*(t_tr/t)*v_0 = (r1/r2)*(t_tr/t)*(v/r4) at temperature t_tr in the transition 
      region, where t is the average coronal temperature.

USAGE:
(set in usage field of opt parameter structure. See input list above.
(1)--include transition region DEM
     additional outputs: dem_tr, dem_cor, logtdem
(2)--exclude transition region DEM (faster)
(3)--include nonthermal electron energy flux
     additional inputs: flux_nt, energy_nt
(4)--compute rad_ratio (25% more computing time)
     additional outputs: dem_tr, dem_cor, logtdem, f_ratio, rad_ratio 

************************************************************************************/

#include "ebtel_functions.h"

int main (void)
{	
	//Use clock to time the entire EBTEL program
	clock_t time_start;
	clock_t time_finish;
	
	//Start the timer
	time_start = clock();
	
	/************************************************************************************
								Variable Declarations 
	************************************************************************************/
	
	/******Variable Declarations******/
	//Struct
	struct ebtel_params_st *params_final;		//Declare instance of structure ebtel_params_st
	struct Option opt;
	
	//Global definitions (declarations in ebtel_functions.h)
	K_B = 1.38e-16;
	KAPPA_0 = 8.12e-7;
	M_P = 1.67e-24;
	PI = 3.14;
	TWO_SEVENTHS = 0.285714;
	SEVEN_HALVES = 3.5;
	TWO_THIRDS = 0.66666667;
	
	//Int
	int i;
	int n;
	int total_time;
	int heating_shape;
	int loop_length;
	
	//Double
	double t_scale;
	double L;
	double h_nano;
	double t_pulse_half;

	FILE *in_file;
	
	//Pointers
	double time_ptr;
	double heat_ptr;
	
	//Char
	char filename_in[64];
	
	//Arrays
	double time[n];
	double heat[n];
	
	/**********************************
	Read in data from parameter file
	**********************************/
	
	//Read in parameters from file
	sprintf(filename_in,"ebtel_parameters.txt");
	in_file = fopen(filename_in,"rt");
	if(in_file == NULL)
	{
		printf("Error! Could not open file.\n");
	}
	
	fscanf(in_file,"%d\n%le\n%d\n%d\%d\n%d\n%d\n%d\n%d\n%d\n%le\n%le\n%le\n%le\n",&total_time,&t_scale,&heating_shape,\
	&loop_length,&opt.usage,&opt.rtv,&opt.dem_old,&opt.dynamic,&opt.solver,&opt.mode,&h_nano,&t_pulse_half,&opt.T0,&opt.n0);
	
	fclose(in_file);
	
	/************************************************************************************
									Initial Parameters
	************************************************************************************/
	
	//Guess maximum array size using the timestep
	n = total_time/t_scale;
	
	//Build the time array
	time_ptr = ebtel_linspace(0,total_time,n);
	for(i=0; i<n; i++)
	{
		time[i] = *(time_ptr + i);
	}
	
	//Define loop half-length and change to appropriate units
	L = 1e8*loop_length;	//convert from Mm to cm
	
	//Set members of the Option opt structure
	opt.heating_shape = heating_shape;
	opt.t_pulse_half = t_pulse_half;
	opt.tau = t_scale;
	opt.h_nano = h_nano;
	if(opt.dynamic == 0)
	{
		opt.classical = 1;
	}
	else
	{
		opt.classical = 0;
	}
	opt.energy_nt = 8.01e-8;	//50 keV in ergs
	
	/************************************************************************************
										Heating
	************************************************************************************/
	//Choose which heating model to use
	//1--triangular pulse (recommended, used in Paper I,II)
	//2--square pulse 
	//3--Gaussian pulse
	
	//Call the heating function and have it return the heating array
	heat_ptr = ebtel_heating(time, t_scale, h_nano, t_pulse_half, n, heating_shape);
	for(i=0;i<n;i++)
	{
		heat[i] = *(heat_ptr + i);
	}
	
	//Set loop length appropriately for Gaussian heating pulse. This heating shape only 
	//appropriate for short loops.
	if(heating_shape==3)
	{
		loop_length=25;
		L=1e+8*loop_length;
	}
	
	/************************************************************************************
									Start the Model
	************************************************************************************/
	
	//Print a header to the screen that gives the user input information
	ebtel_print_header(n, heating_shape, loop_length, total_time, opt);
	
	//Make the call to the ebtel_loop_solver function. This function sets the members of the structure params_final. Each member 
	//is a pointer to an array
	params_final = ebtel_loop_solver(n, L, total_time, time, heat, opt);
	
	//Save the heating and time arrays to the parameter structure
	params_final->time = time_ptr;
	params_final->heat = heat_ptr;
	
	/************************************************************************************
									Save the Data
	************************************************************************************/
	
	//Write the contents of params_final to a file. See output for filename.
	ebtel_file_writer(loop_length, n, opt, params_final);
	
	/****************Done writing data to file. Free up memory reserved by pointers.******************/
	
	//Free up memory used by the structure params_final
	ebtel_free_mem(params_final);
	free(heat_ptr);
	free(time_ptr);
	
	//Stop the timer
	time_finish = clock();
	
	//Time elapsed
	printf("The process took %ld seconds to run\n",(time_finish - time_start)/CLOCKS_PER_SEC);
	
	//Exit with no errors
	return 0;
	
}