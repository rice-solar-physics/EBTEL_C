/************************************************************************************

FILENAME: ebtel_main.c

AUTHOR Will Barnes

DATE: created: 7 March 2014

-------------------------------------------------------------------------------------
EBTEL (Enthalpy Based Thermal Evolution of Loops) computes 0-D
hydrodynamic equations.  This software was originally developed by Klimchuk et al,. (2008) and
later improved upon by Cargill et al., (2012a,b). It solves the set of hydrostatic equations by 
averaging over the loop half-length and then integrating in time. Additional details are available
in the references given above.

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
(set in usage field of opt parameter structure. See input list above.)
dem--include transition region DEM
     additional outputs: dem_tr, dem_cor, logtdem (recommended)
no_dem--exclude transition region DEM (faster)
nt_ebeam--include nonthermal electron energy flux
     additional inputs: flux_nt, energy_nt
rad_ratio--compute rad_ratio (requires longer compute time)
     additional outputs: dem_tr, dem_cor, logtdem, f_ratio, rad_ratio 

************************************************************************************/

#include "ebtel_functions.h"

int main (int argc, char *argv[])
{	
	//Use clock to time the entire EBTEL program
	clock_t time_start;
	clock_t time_diff;
	double time_elapsed;
	
	//Start the timer
	time_start = clock();
	
	/************************************************************************************
								Variable Declarations 
	************************************************************************************/
	
	/******Variable Declarations******/
	//Struct
	struct ebtel_params_st *params_final;
	struct Option *opt;							
	
	//Global definitions (declarations in ebtel_functions.h)
	//KAPPA_0 = 1e-6;
	KAPPA_0 = 8.12e-7;
	PI = 3.14159265359;
	TWO_SEVENTHS = 2./7.;
	SEVEN_HALVES = 3.5;
	TWO_THIRDS = 2./3.;
	M_EL = 9.11e-28;
	
	//Set global variables based on He/H abundance
	ebtel_calc_abundance();
	
	//Int
	int i,n;
	int quiet_flag = 0;
	double L;
	char filename_in[250];
	double *kptr;
	
	/**********************************
	Read configuration file
	**********************************/

	//Read in parameters from file using xmllib library
	//Set default filename
	sprintf(filename_in,"../config/ebtel_config.xml");
	//Check if a filename was specified at command line
	for(i = 0; i<argc; i++)
	{
		//Read in filename
		if(strcmp(argv[i],"quiet")!=0 && i>0)
		{
			sprintf(filename_in,"%s",argv[i]);
		}
		else if(strcmp(argv[i],"quiet")==0)
		{
			//Raise the quiet flag
			quiet_flag = 1;
		}
	}
	//Pass filename to struct setter to read inputs into opt structure
	opt = ebtel_input_setter(filename_in);
	
	
	/************************************************************************************
									Initial Parameters
	************************************************************************************/
	
	//Set total number of steps using the initial timestep and total time
	//When using the adaptive method, this can be increased to avoid segmentation fault runtime error.
	//Set total number of steps using the initial timestep and total time
	if(strcmp(opt->solver,"euler")==0 || strcmp(opt->solver,"rk4")==0)
	{
		//For static timesteps, this is just the total time divided by the timestep
		n = ceil(opt->total_time/opt->tau)+1;
	}
	else if(strcmp(opt->solver,"rka4")==0)
	{
		//When using the adaptive method, this is a guess and additional memory will be allocated if necessary
 	   	n = ceil(opt->total_time);
	}
	else
	{
		printf("Invalid solver option.\n");
		exit(0);
	}
	
	//Define loop half-length and change to appropriate units
	L = 1e8*opt->loop_length;	//convert from Mm to cm
	
	//Set non-thermal electron energy structure option
	opt->energy_nt = 8.01e-8;	//50 keV in ergs
	
	//Set temperature bins for calculating radiative loss function
 	//Make the kpar array
	NK=6;	//There is a correspondence with the length of KPAR[]; if NK ever changes, change length of KPAR[] in header file
	kptr = ebtel_kpar_set(opt->rad_option);
 	for(i=0; i<NK; i++)
 	{
 		KPAR[i] = *(kptr + i);
 	}
	free(kptr);
	kptr=NULL;
	
	/************************************************************************************
									Heating
	************************************************************************************/
	
	//Configure start times, end times, and amplitudes of heating events either from input
	//from above parameters, through normally distributed start times and amplitudes following
	//a power-law distribution or by reading in all three parameters from specified input
	//files.
	ebtel_heating_config(opt,filename_in);
	
	/************************************************************************************
									Start the Model
	************************************************************************************/
	
	//Print a header to the screen that gives the user input information
	//(If you're doing a large parameter sweep, use the 'quiet' command line option.)
	if(quiet_flag == 0)
	{
		ebtel_print_header(n, opt);
	}
	
	//Make the call to the ebtel_loop_solver function. This function sets the members of the structure params_final. Each member 
	//is a pointer to an array
	params_final = ebtel_loop_solver(n, L, opt);
	
	/************************************************************************************
									Save the Data
	************************************************************************************/
	
	//Write the contents of params_final to a file. See output for filename.
	ebtel_file_writer(opt, params_final);
	
	//Count the number of events and print it to the screen
	int num_q_events = ebtel_count_events(params_final,opt);
	printf("Number of simulated heating events: %d\n",num_q_events);
	
	/****************Done writing data to file. Free up memory reserved by pointers.******************/
	
	//Free up memory used by the structures params_final and opt 
	ebtel_free_mem(params_final,opt);
	
	//Stop the timer
	time_diff = clock() - time_start;
	time_elapsed = time_diff*1000/CLOCKS_PER_SEC;
	
	//Time elapsed
	printf("The process took %f milliseconds to run.\n",time_elapsed);
	
	//Exit with no errors 
	return 0;
	
}