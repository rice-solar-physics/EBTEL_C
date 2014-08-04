/***********************************************************************************

FILENAME: ebtel_functions_heating.c

AUTHOR Will Barnes

DATE: created: 4 August 2014

DESCRIPTION: This file contains all of the functions relating to the ad-hoc heating
imposed in the EBTEL model. 

***********************************************************************************/

//Include appropriate header file
#include "ebtel_functions.h"

/***********************************************************************************

FUNCTION NAME: ebtel_heating

FUNCTION_DESCRIPTION: This function sets up the heating array to be used in the heatings
of our coronal loop. Currently, three different types of heating are available: triangular,
square, or Gaussian pulse. 

INPUTS:
	t--time array for our model
	opt--structure that provides all necessary input parameters
	
OUTPUTS:
	heat--heating at time time

***********************************************************************************/

double ebtel_heating(double t, struct Option *opt)
{
	//Declare variables
	int i;
	double h_back;
	double h_thick;
	double t_pulse;
	double heat;

	//First set some general parameters
	h_back = 3.4e-6;
	h_thick = 0;
	t_pulse = 2*opt->t_pulse_half;
		
	//Set the heating as the background heating
	heat = h_back;
	
	//Check all heating intervals to see if we have fallen into one of them
	//Test all timing intervals
	for(i=0;i<opt->num_events;i++)
	{
		//Check if we are inside the heating pulse interval
		if(t >= *(opt->t_start_array + i) && t <= (*(opt->t_start_array + i) + t_pulse) )
		{
			//If so, call the heating profile function to generate the correct pulse
			heat = ebtel_heating_profile(t,*(opt->t_start_array + i),*(opt->amp + i),opt);
			heat = heat + h_back;
		}
	}
	
	//Return the heating value
	return heat;
}

/***********************************************************************************

FUNCTION NAME: ebtel_heating_profiles

FUNCTION_DESCRIPTION: This function chooses the heating profile from the heating_shape
input and returns a value for the heating based on the selected profile and the current
time.

INPUTS:
	t--time array for our model
	t_start--starting time of the current heating pulse
	h_nano--amplitude of the current heating pulse
	opt--structure that provides all necessary input parameters
	
OUTPUTS:
	heat--heating at time time

***********************************************************************************/

double ebtel_heating_profile(double t, double t_start, double h_nano, struct Option *opt)
{
	//Variable declarations and definitions
	double t_pulse = 2*opt->t_pulse_half;
	double t_mid = t_start + t_pulse/2.;
	double t_end = t_start + t_pulse;
	double heat;
	
	//Choose which heating model to use
	//1--triangular pulse (recommended, used in Paper I,II)
	//2--square pulse 
	//3--Gaussian pulse
	//Additional heating functions should be added here.
	
	if(opt->heating_shape == 1)
	{
		//Triangular Pulse
		if(t < t_mid)
		{
			heat = h_nano*(t - t_start)/(t_pulse/2.);
		}
		else 
		{
			heat = -h_nano*(t - t_end)/(t_pulse/2.);
		}

    }
	else if(opt->heating_shape == 2)
	{
		//Square pulse
		heat = h_nano;
		
    }
	else if(opt->heating_shape == 3)
	{
		//Gaussian
		heat = h_nano*exp(-pow((t - t_mid),2)/(2*pow(opt->t_pulse_half,2)));
	}
	else
	{
		printf("Invalid heating profile choice. Exiting program.\n");
		exit(0);
	}
	
	//Return the resulting heating value
	return heat;
}

/***********************************************************************************

FUNCTION NAME: ebtel_rand_limit

FUNCTION_DESCRIPTION: This function gives a uniformly distributed random number from
zero to limit where limit is specified by the input argument. 

INPUTS:
	limit--maximum of uniform distribution
	
OUTPUTS:
	retval--resulting random number

***********************************************************************************/

double ebtel_rand_limit(double limit)
{
	double divisor = RAND_MAX/limit;
	double retval;
	
	do{
		retval = rand()/divisor;
	}while(retval > limit);
	
	return retval;
}

/***********************************************************************************

FUNCTION NAME: ebtel_box_muller

FUNCTION_DESCRIPTION: This function is an implementation of the Box-Muller algorithm
to achieve a normal distribution from a random distribution. A description of the 
algorithm can be found in 'Numerical Recipes' by Press though this implementation is 
not equivalent. The function returns a structure holding two resulting values as well
as a flag telling whether an old value is available for use.

INPUTS:
	x1--uniformly distributed random number
	x2--uniformly distributed random number
	save--leftover value
	flag--flag to indicate leftovers
	
OUTPUTS:
	bm_st--structure holding normally distributed values and updated flag

***********************************************************************************/

struct box_muller_st *ebtel_box_muller(double x1,double x2,double save,int flag)
{
	//Declare variables
	struct box_muller_st *bm_st = malloc(sizeof(struct box_muller_st));
	double z1,z2;
	
	//Check if the flag has been raised. If not, compute two values
	if(flag==1)
	{
		bm_st->z = save;
		bm_st->z_save = 0.;
		bm_st->flag = 0;
	}
	else
	{
		z1 = sqrt(-2.*log(x1))*cos(2.*PI*x2);
		z2 = sqrt(-2.*log(x1))*sin(2.*PI*x2);
		bm_st->z = z1;
		bm_st->z_save = z2;
		bm_st->flag = 1;
	}
	
	//Return the structure
	return bm_st;
}

/***********************************************************************************

FUNCTION NAME: ebtel_power_law

FUNCTION_DESCRIPTION: This function takes in uniformly distributed random number
and produces a power-law distribution over range [x0,x1] with index alpha.

INPUTS:
	x0--lower limit to power-law distribution
	x1--upper limit to power-law distribution
	y--uniformly distributed random number
	alpha--power-law index
	
OUTPUTS:
	x--random number that follows power law distribution

***********************************************************************************/

double ebtel_power_law(double x0, double x1, double y, double alpha)
{
	//Declare variables
	double x;
	double term1,term2;
	
	//Calculate first term
	term1 = pow(x1,(alpha + 1.)) - pow(x0,(alpha + 1.));
	term2 = pow(x0,(alpha + 1.));
	x = pow((term1*y + term2),(1/(alpha + 1)));
	
	//Return the result
	return x;
}

/***********************************************************************************

FUNCTION NAME: ebtel_bubble_sort

FUNCTION_DESCRIPTION: This function is an implementation of the well-known bubble-
sort algorithm and sorts the array in ascending order.

INPUTS:
	array--array of numbers to be sorted
	array_length--length of the array to be sorted
	
OUTPUTS:
	sorted_array--pointer to sorted array
***********************************************************************************/

double * ebtel_bubble_sort(double array[],int array_length)
{
	//Declare variables
	int switch_flag;
	int i;
	double tmp;
	double *sorted_array = malloc(sizeof(double[array_length]));
	
	//Begin first loop that won't end until sorting is done
	do{
		//Reset the flag
		switch_flag = 0;
		
		//Loop over array
		for(i=0; i<(array_length-1); i++)
		{
			if(array[i] > array[i+1])
			{
				//Swap the values
				tmp = array[i+1];
				array[i+1] = array[i];
				array[i] = tmp;
				
				//Raise the flag
				switch_flag = 1;
				
				//Break out of the inner loop
				break;
			}
		}
		
	}while(switch_flag == 1);
	
	//Make a pointer to the new array
	for(i = 0; i<array_length; i++)
	{
		sorted_array[i] = array[i];
	}
	
	//Return sorted array
	return sorted_array;
}
