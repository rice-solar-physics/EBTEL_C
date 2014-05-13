/***********************************************************************************

FILENAME: ebtel_functions_util.c

AUTHOR Will Barnes

DATE: created: 7 April 2014

DESCRIPTION: This file contains utility functions needed in EBTEL. For additional 
details, see ebtel_main.c.

***********************************************************************************/

//Include appropriate header file
#include "ebtel_functions.h"

/***********************************************************************************

FUNCTION NAME: ebtel_print_header

FUNCTION_DESCRIPTION: This function prints a header to the screen and tells the user
what input is being used in this implementation of the model. 

INPUTS:
	a--starting point
	b--ending point
	n--number of entries
	
OUTPUTS:
	

***********************************************************************************/

void ebtel_print_header(int n, int heating_shape, int loop_length, int total_time, struct Option opt)
{
	//Print a header and tell the user what options are being used to begin the model
	printf("************************************************************************************\n");
	printf("           Enthalpy Based Thermal Evolution of Loops (EBTEL)                        \n");
	printf("************************************************************************************\n\n");
	printf("Original code written in IDL by J.A. Klimchuk, S. Patsourakos, P.J. Cargill\n");
	printf("See Klimchuk, J.A, S. Patsourakos & P.J. Cargill 2008, ApJ 682:1351-2362\n");
	printf("See also Cargill, P.J., S.J. Bradshaw & J.A. Klimchuk 2012, ApJ 752:161-174\n\n");
	printf("Translation into the C Programming Language by Will Barnes,\nDept. of Physics & Astronomy, Rice University (2014)\n");
	printf("************************************************************************************\n\n");
	printf("INPUTS\n");
	printf("------\n");
	printf("Number of steps: %d\n",n);
	printf("Total time: %d s\n",total_time);
	printf("Loop half-length: %d Mm\n",loop_length);
	printf("Usage option(see documentation): %d\n",opt.usage);
	if(heating_shape==1)
	{printf("Heating: Triangular heating pulse\n");
	}
	else if(heating_shape==2)
	{printf("Heating: Square heating pulse\n");
	}
	else
	{printf("Heating: Gaussian heating pulse (note: L has been altered for this heating option)\n");
	}
	if(opt.solver==1)
	{printf("Solving equations using fourth order Runge-Kutta routine\n");
	}
	else if(opt.solver==2)
	{printf("Solving equations using adaptive fourth order Runge-Kutta routine\n");
	}
	else 
	{printf("Solving equations using Euler method\n");
	}
	if(opt.rtv==1)
	{printf("Using Rosner-Tucker-Vaiana Loss Function\n");
	}
	else
	{printf("Using Raymond-Klimchuk Loss Function\n");
	}
	if(opt.dynamic==1)
	{printf("Using dynamic heat flux calculation\n");
	}
	else
	{printf("Using classical heat flux calculation\n");
	}
	if(opt.usage==1 || opt.usage==4)
	{
		if(opt.dem_old==1)
		{printf("Using old method to calculate DEM in the TR\n");
		}
		else
		{printf("Using new method to calculate DEM in the TR\n");
		}
	}
	printf("\n");
}

/***********************************************************************************

FUNCTION NAME: ebtel_file_writer

FUNCTION_DESCRIPTION: This function prints the data from the structure ebtel_params_st
to a file. It produces separate files for the loop data and DEM data. 

INPUTS:
	loop_length--loop half length in Mm
	n--length of the arrays to be printed
	opt--Option structure with inputs necessary for labeling
	params_final--pointer to structure of output data
	
OUTPUTS:
	none

***********************************************************************************/

void ebtel_file_writer(int loop_length, int n, struct Option opt, struct ebtel_params_st *params_final)
{
	//Declare variables
	int i;
	double f_ratio[n];
	double rad_ratio[n];
	char filename_out[64];
	char filename_out_dem[64];
	FILE *out_file;
	
	//Open the file that we are going to write the data to 
	sprintf(filename_out,"data/ebteldatL%du%dh%ds%d.txt",loop_length,opt.usage,opt.heating_shape,opt.solver);
	out_file = fopen(filename_out,"wt");
	
	//Tell the user where the results were printed
	printf("The results were printed to the file %s\n",filename_out);
	
	//The members of the structure params_final have now been set so we need to unpack them and set our arrays so that we can
	//easily save our data.
	for(i = 0; i<n; i++)
	{	
		//If we used usage = 4, then we need to save f_ratio and rad_ratio as well.
		//We set them to zero otherwise just as a placeholder.
		if(opt.usage == 4)
		{
			rad_ratio[i] = *(params_final->rad_ratio + i);
			f_ratio[i] = *(params_final->f_ratio + i);
		}
		else
		{
			rad_ratio[i] = 0;
			f_ratio[i] = 0;
		}
		
		//Print the data to the file filename using tab delimited entries
		fprintf(out_file,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",*(params_final->time + i),*(params_final->temp + i),*(params_final->ndens + i),*(params_final->press + i),*(params_final->vel + i),*(params_final->tapex + i),*(params_final->napex +i),*(params_final->papex + i),*(params_final->cond + i),*(params_final->rad_cor + i),rad_ratio[i],f_ratio[i],*(params_final->heat + i),*(params_final->coeff_1 + i));
		
	}
	
	//Close the file when we are done writing data to it
	fclose(out_file);
	
	//If we chose to calculate the TR DEM, we need to write this data to a separate file.
	if(opt.usage==1 || opt.usage==4)
	{
		//Make the DEM data filename
		sprintf(filename_out_dem,"data/ebteldemdatL%du%dh%ds%d.txt",loop_length,opt.usage,opt.heating_shape,opt.solver);
		
		//Tell the user where the DEM data was printed to
		printf("The DEM results were printed to the file %s\n",filename_out_dem);
		
		//Open the DEM data file
		out_file = fopen(filename_out_dem,"wt");
		
		for(i=0; i<451; i++)
		{	
			//Now write this data to a file. This will be a secondary table underneath the first
			fprintf(out_file,"%e\t%e\t%e\t%e\n",*(params_final->logtdem + i), *(params_final->dem_tr_log10mean + i), \
			*(params_final->dem_cor_log10mean + i),*(params_final->dem_tot_log10mean + i));
		}
		
		//Close the DEM data file
		fclose(out_file);
	}
	
}

/***********************************************************************************

FUNCTION NAME: ebtel_linspace

FUNCTION_DESCRIPTION: This function is similar to the linspace function in MATLAB. It creates a one-dimensional array with specified number of entries and starting and ending point. 

INPUTS:
	a--starting point
	b--ending point
	n--number of entries
	
OUTPUTS:
	linspace--pointer to array with n points starting at a and ending at b

***********************************************************************************/

double * ebtel_linspace( int a, int b, int n)
{
	int i;
	double *linspace = malloc(sizeof(double[n]));
	
	//Declare necessary variables
	double interval = (b - a)/(n-1);	//spacing between points
	
	//Make the array
	for(i = 0; i<n; i++)
	{
		linspace[i] = a + i*interval;
	}

	return linspace;
}

 /**********************************************************************************
 
 Function name: ebtel_avg_val
 
 Function description: This function accepts an array and calculates its mean or average
 value. 
 
 Input
	double numbers []	Array of double values
	int length			Integer value, length of array numbers
 
 Return
	double mean			Double value returned by the function.
 
 *********************************************************************************/
 
 double ebtel_avg_val (double numbers[], int length)
 {
 	//Declare some variables
 	int i;
 	double mean;
 	double sum = 0;
 	
 	//Calculate the sum
 	for (i = 0; i<length; i++)
 	{
 		sum = sum + numbers[i];	
 	}
 	
 	//Calculate the mean
 	mean = sum/length;
 	
 	//Return the mean
 	return mean;
 }
 
/**********************************************************************************
 
 Function name: ebtel_max_val
 
 Function description: This function accepts two numbers and returns the maximum of 
 those two numbers.
 
 Input
	double num_1		first number to compare
	double num_2		second number to compare	
 
 Return
	double max_val		maximum of num_1 and num_2
 
 *********************************************************************************/
 
 double ebtel_max_val(double num_1, double num_2)
 {
 	//Declare max value
 	double max_val;
 	
 	if(num_1 > num_2)	
 	{
 		max_val = num_1;
 	}
 	else if(num_2 > num_1)
 	{
 		max_val = num_2;
 	}
 	else					//Consider case where num_1 = num_2; Could set either one to the max since it does not matter
 	{
 		max_val = num_1;
 	}
 	
 	return max_val;
 }
 
 /**********************************************************************************
 
 Function name: ebtel_min_val
 
 Function description: This function accepts two numbers and returns the minimum of
 these two numbers.
 
 Input
	double num_1		first number to compare
	double num_2		second number to compare	
 
 Return
	double min_val		minimum of num_1 and num_2
 
 *********************************************************************************/
 
 double ebtel_min_val(double num_1, double num_2)
 {
 	//Declare min value
 	double min_val;
 	
 	if(num_1 < num_2)
 	{
 		min_val = num_1;
 	}
 	else if(num_2 < num_1)
 	{
 		min_val = num_2;
 	}
 	else					//Consider the case where the two are equal. Doesn't matter which we return
 	{
 		min_val = num_1;
 	}
 	
 	return min_val;
 }
 
 /**********************************************************************************
 
 Function name: ebtel_free_mem
 
 Function description: This function frees up memory that was malloc'd by the ebtel_loop
 _solver function in order to return the param structure
 
 Input
	struct ebtel_params_st *par_struct		pointer to structure memory that will be freed	
 
 Return
	none
 
 *********************************************************************************/
 void ebtel_free_mem(struct ebtel_params_st *par_struct)
 {
 	//First check that the pointer is valid
 	assert(par_struct != NULL);
 	
 	//Free the memory of each of the structure members as well as the structure itself
 	free(par_struct->time);
 	par_struct->time = NULL;
 	free(par_struct->heat);
 	par_struct->heat = NULL;
 	free(par_struct->temp);
 	par_struct->temp = NULL;
	free(par_struct->ndens);
	par_struct->ndens = NULL;
	free(par_struct->press);
	par_struct->press = NULL;
	free(par_struct->vel);
	par_struct->vel = NULL;
	free(par_struct->tapex);
	par_struct->tapex = NULL;
	free(par_struct->napex);
	par_struct->napex = NULL;
	free(par_struct->papex);
	par_struct->papex = NULL;
	free(par_struct->coeff_1);
	par_struct->coeff_1 = NULL;
	free(par_struct->logtdem);
	par_struct->logtdem = NULL;
	free(par_struct->f_ratio);
	par_struct->f_ratio = NULL;
	free(par_struct->rad_ratio);
	par_struct->rad_ratio = NULL;
	free(par_struct->cond);
	par_struct->cond = NULL;
	free(par_struct->rad_cor);
	par_struct->rad_cor = NULL;
	free(par_struct->dem_tr_log10mean);
	par_struct->dem_tr_log10mean = NULL;
	free(par_struct->dem_cor_log10mean);
	par_struct->dem_cor_log10mean = NULL; 
	free(par_struct->dem_tot_log10mean);
	par_struct->dem_tot_log10mean = NULL;
	
	//Free memory reserved for the structure
	free(par_struct);
 }
 