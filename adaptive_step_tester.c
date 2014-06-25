//adaptive_step_tester.c
//Will Barnes
//19 June 2014

//This function will test the adaptive method used in the EBTEL model by feeding it various test functions and then outputing the timestep.

#include "ebtel_functions.h"

//NOTE: Put all of the rk derivs junk in the rk derivs file to eliminate error. Then call the adaptive rk method from main

//Derivs routine for cubic function
double * rk_derivs_cubic(double state[],double time)
{
	//Declare variables
	double x;
	double dxdt;
	double *derivs = malloc(sizeof(double[1]));
	
	//Unravel the state vector
	x = state[0];
	
	//Compute the derivative
	dxdt = 1.;
	
	derivs[0] = dxdt;
	
	return derivs;
}

//Driver routine for testing 
int main(void)
{
	//Decalarations for printing to a file
	FILE *out_file;
	
	//Declare variables
	double tau = .1;
	double total_time = 1000.;
	double start_time = -1000.;
	int i = 0;
	int j;
	int ar_size = (total_time-start_time)/tau;
	
	//Declare arrays used for saving
	double pos[ar_size];
	double time[ar_size];
	double time_step[ar_size];
	double temp[1];
	double space_step;
	
	//Declare instances of par and opt structures to appease rka function
	struct rk_params par;
	struct Option opt;
	
	//Declare structure to return adapt information
	struct ebtel_rka_st *adapt;
	
	//Initialize arrays
	pos[0] = pow(start_time,2.);
	time[0] = start_time;
	time_step[0] = tau;
	temp[0] = pos[0];
	
	//Check and see if directory 'data' exists. If it does not, then create a new one.
	struct stat st = {0};
	if(stat("test_data",&st) == -1){
		mkdir("test_data",0777);
	}
	
	//Open the file
	out_file = fopen("test_data/adaptive_testing.txt","wt");
	
	//Begin loop over time
	while(time[i]<total_time)
	{
		
		//Call the adaptive method
		adapt = ebtel_rk_adapt(temp,1,time[i],tau,1e-10,par,opt);
		tau = adapt->tau;
		
		//Update vectors
		pos[i+1] = *(adapt->state + 0);
		time[i+1] = time[i] + tau;
		time_step[i+1] = tau;
		temp[0] = pos[i+1];
					
		//Free the structure as it will be malloc'd on the next go around
		free(adapt->state);
		adapt->state = NULL;
		free(adapt);
		
		
		/*
		//Use an Euler method
		time[i+1] = time[i] + tau;
		temp = rk_derivs_cubic()
		space_step = tau*2*pow(time[i],1.);
		//Print position before stepping
		pos[i+1] = pos[i] + space_step;
		time_step[i+1] = tau;
		*/
		
		printf("pos[%d] = %f\n",i,pos[i]);
		fprintf(out_file,"%f\t%f\t%f\n",time[i],pos[i],time_step[i]);
		
		//Increment counter
		i++;
	}
	
	//Close the file
	fclose(out_file);
	
	return 0;
	
}