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
	clock_t time_diff;
	double time_elapsed;
	
	//Start the timer
	time_start = clock();
	
	/************************************************************************************
								Variable Declarations 
	************************************************************************************/
	
	/******Variable Declarations******/
	//Struct
	struct ebtel_params_st *params_final;		//Declare instance of structure ebtel_params_st
	struct Option *opt = malloc(sizeof(struct Option));
	struct box_muller_st *bm_st;
	
	//Global definitions (declarations in ebtel_functions.h)
	//KAPPA_0 = 1e-6;
	KAPPA_0 = 8.12e-7;
	PI = 3.14159265359;
	TWO_SEVENTHS = 2./7.;
	SEVEN_HALVES = 3.5;
	TWO_THIRDS = 2./3.;
	
	//Set global variables based on He/H abundance
	ebtel_calc_abundance();
	
	//Int
	int n;
	int i;
	int heating_shape;
	int loop_length;
	int num_events;
	int alpha;
	int bm_flag = 0;
	
	//Double
	double total_time;
	double t_scale;
	double L;
	double h_nano;
	double t_pulse_half;
	double t_start;
	double mean_t_start,std_t_start;
	double amp_0,amp_1;
	double x1,x2;
	double tmp,save;
	double limit = 1.;
	
	FILE *in_file;
	FILE *in_file_start;
	FILE *in_file_amp;
	FILE *in_file_end;
	double *sort_ptr1;
	double *sort_ptr2;
	
	//Char
	char filename_in[64];
	char t_start_switch[64];
	char amp_switch[64];
	char t_end_switch[64];
	char start_file[64];
	char end_file[64];
	char amp_file[64];
	
	/**********************************
	Read in data from parameter file
	**********************************/
	
	//Read in parameters from file
	sprintf(filename_in,"ebtel_parameters.txt");
	in_file = fopen(filename_in,"rt");
	if(in_file == NULL)
	{
		printf("Error! Could not open file.\n");
		return 1;
	}
	
	fscanf(in_file,"%le\n%le\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%le\n%le\n%le\n%d\n%le\n%le%le\n",&total_time,&t_scale,&heating_shape,&loop_length,&opt->usage,&opt->rtv,&opt->dem_old,&opt->dynamic,&opt->solver,&opt->mode,&h_nano,&t_pulse_half,&t_start,&opt->index_dem,&opt->error,&opt->T0,&opt->n0);
	
	fclose(in_file);
	
	/************************************************************************************
									Initial Parameters
	************************************************************************************/
	
	//Set total number of steps using the initial timestep and total time
	//When using the adaptive method, this can be increased to avoid segmentation fault runtime error.
	n = ceil(2*total_time/t_scale);
	
	//Define loop half-length and change to appropriate units
	L = 1e8*loop_length;	//convert from Mm to cm
	
	//Set members of the Option opt structure
	opt->heating_shape = heating_shape;
	opt->t_pulse_half = t_pulse_half;
	opt->t_start = t_start;
	opt->tau = t_scale;
	opt->h_nano = h_nano;
	opt->energy_nt = 8.01e-8;	//50 keV in ergs
	
	/************************************************************************************
									Heating
	************************************************************************************/

	//Calculate start times and amplitudes from appropriate distributions if we are using 

	//Read in input parameters from heating input file
	sprintf(filename_in,"ebtel_heating_parameters.txt");
	in_file = fopen(filename_in,"rt");
	if(in_file == NULL)
	{
		printf("Error! Could not open heating parameters file.\n");
		return 1;
	}
	fscanf(in_file,"%d\n%le\n%le\n%d\n%le\n%le\n%s\n%s\n%s\n%s\n%s\n%s\n",&num_events,&mean_t_start,&std_t_start,&alpha,&amp_0,&amp_1,t_start_switch,amp_switch,t_end_switch,start_file,amp_file,end_file);
	fclose(in_file);

	//Set the number of heating events in the input structure
	opt->num_events = num_events;

	//Declare amplitude and start time arrays
	double amp[num_events];
	double t_start_array[num_events];
	double t_end_array[num_events];

	//Reserve memory for amplitude and start time arrays in opt structure
	opt->t_start_array = malloc(sizeof(double[num_events]));
	opt->amp = malloc(sizeof(double[num_events]));
	opt->t_end_array = malloc(sizeof(double[num_events]));

	//Seed the random number generator
	srand(time(NULL));

	//Calculate the start times and amplitudes
	//Begin loop to set start times
	for(i=0;i<num_events;i++)
	{
		//Set random numbers for either start time or amplitudes
		if(strcmp(t_start_switch,"random") ==0 || strcmp(amp_switch,"random") == 0)
		{
			//Initialize the two random variables
			x1 = ebtel_rand_limit(limit);
			x2 = ebtel_rand_limit(limit);
		}
	
		//Use uniformly or normally distributed start times
		if(strcmp(t_start_switch,"uniform") == 0)
		{
			//Start times separated by two pulse durations (following Reep et al. 2013)
			t_start_array[i] = t_start + 2.*i*(2*t_pulse_half);
		}
		else if(strcmp(t_start_switch,"random") == 0)
		{
			//Use the Box-Muller method to do the normal distribution
			bm_st = ebtel_box_muller(x1,x2,save,bm_flag);
			tmp = bm_st->z;
			save = bm_st->z_save;
			bm_flag = bm_st->flag;

			//Save the 'denormalized' normally distributed start time
			t_start_array[i] = std_t_start*tmp + mean_t_start;
		
			//Free the structure
			free(bm_st);
			bm_st = NULL;
		}
		else if(strcmp(t_start_switch,"file") == 0)
		{
			//Open file on the first iteration
			if(i==0)
			{
				//DEBUG--print the file name
				printf("Heating start time file: %s\n",start_file);
				in_file_start = fopen(start_file,"rt");
				if(in_file_start==NULL)
				{
					printf("Error! Could not open heating start time file.\n");
					return 1;
				}
			}
			
			//Read in start times from file
			fscanf(in_file_start,"%le\n",&t_start_array[i]);
				
			//Close file on last iteration 
			if(i==(num_events-1))
			{
				fclose(in_file_start);
			}
		}
		else
		{
			printf("Invalid heating start time option. Choose either uniform, file or normally distributed\n");
			exit(0);
		}
	
		//Use uniform amplitudes, amplitudes given by power law distribution, or read them in from a file
		if(strcmp(amp_switch,"uniform") == 0)
		{
			amp[i] = h_nano;
		}
		else if(strcmp(amp_switch,"random") == 0)
		{
			//Compute the amplitude according to a power-law distribution
			amp[i] = ebtel_power_law(amp_0,amp_1,x1,alpha);
		}
		else if(strcmp(amp_switch,"file") == 0)
		{
			//Open file on the first iteration
			if(i==0)
			{
				in_file_amp = fopen(amp_file,"rt");
				if(in_file_amp==NULL)
				{
					printf("Error! Could not open heating amplitude file.\n");
					return 1;
				}
			}
			
			//Read in start times from file
			fscanf(in_file_amp,"%le\n",&amp[i]);
				
			//Close file on last iteration 
			if(i==(num_events-1))
			{
				fclose(in_file_amp);
			}
		}
		else
		{
			printf("Invalid heating amplitude option. Choose either uniform, file or power-law distribution\n");
			exit(0);
		}
		
		//Use uniform pulse times (as defined by configuration file) or read in pulse times from separate file
		if(strcmp(t_end_switch,"uniform") == 0)
		{
			//Set array of pulse times from configuration file
			t_end_array[i] = 2*t_pulse_half + t_start_array[i];
		}
		else if(strcmp(t_end_switch,"file") == 0)
		{
			//Open file on the first iteration
			if(i==0)
			{
				in_file_end = fopen(end_file,"rt");
				if(in_file_end==NULL)
				{
					printf("Error! Could not open heating end time file.\n");
					return 1;
				}
			}
			
			//Read in start times from file
			fscanf(in_file_end,"%le\n",&t_end_array[i]);
				
			//Close file on last iteration 
			if(i==(num_events-1))
			{
				fclose(in_file_end);
			}
		}
		else
		{
			printf("Invalid heating pulse option. Choose either uniform or file option\n");
			exit(0);
		}

	}

	//If the start times are random, Sort start and end times in ascending order and set pointers in opt structure
	if(strcmp(t_start_switch,"random") == 0)
	{
		sort_ptr1 = ebtel_bubble_sort(t_start_array,num_events);
		sort_ptr2 = ebtel_bubble_sort(t_end_array,num_events);	
		for(i=0;i<num_events;i++)
		{
			t_start_array[i] = *(sort_ptr1 + i);
			t_end_array[i] = *(sort_ptr2 + i);			
		}
		free(sort_ptr1);
		sort_ptr1=NULL;
		free(sort_ptr2);
		sort_ptr2=NULL;
	}
	
	//Save the start time, amplitude, and pulse arrays to the opt structure
	for(i=0; i<num_events; i++)
	{
		opt->t_start_array[i] = t_start_array[i];
		opt->t_end_array[i] = t_end_array[i];
		opt->amp[i] = amp[i];
	}

	
	/************************************************************************************
									Start the Model
	************************************************************************************/
	
	//Print a header to the screen that gives the user input information
	//(If you're doing a large parameter sweep, this line should be commented out.)
	ebtel_print_header(n, heating_shape, loop_length, total_time, opt);
	
	//Make the call to the ebtel_loop_solver function. This function sets the members of the structure params_final. Each member 
	//is a pointer to an array
	params_final = ebtel_loop_solver(n, L, total_time, opt);
	
	/************************************************************************************
									Save the Data
	************************************************************************************/
	
	//Write the contents of params_final to a file. See output for filename.
	ebtel_file_writer(loop_length, opt, params_final);
	
	/****************Done writing data to file. Free up memory reserved by pointers.******************/
	
	//Free up memory used by the structure params_final
	ebtel_free_mem(params_final);
	//Free the t_start and amp arrays
	free(opt->t_start_array);
	free(opt->amp);
	free(opt->t_end_array);
	opt->t_start_array = NULL;
	opt->amp = NULL;
	opt->t_end_array = NULL;
	//Free the memory used by opt input structure
	free(opt);
	
	//Stop the timer
	time_diff = clock() - time_start;
	time_elapsed = time_diff*1000/CLOCKS_PER_SEC;
	
	//Time elapsed
	printf("The process took %f milliseconds to run\n",time_elapsed);
	
	//Exit with no errors 
	return 0;
	
}