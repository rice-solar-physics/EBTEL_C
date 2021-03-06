/***********************************************************************************

FILENAME: ebtel_functions_solvers.c

AUTHOR Will Barnes

DATE: created: 19 April 2014

***********************************************************************************/

//Include appropriate header file
#include "ebtel_functions.h"

/**********************************************************************************

FUNCTION NAME: ebtel_derivs

FUNCTION DESCRIPTION: This function computes the needed derivatives for both the 
Euler and Runge-Kutta solvers in EBTEL. In the Euler case, the resulting updated 
parameters are passed back to the main loops function. For the RK4 cases (adaptive 
and non-adaptive), the resulting derivatives are passed to the relevant RK functions.
 
INPUTS:
	s--vector that stores the current state of the system
	t--current time
	par--structure that holds necessary parameters
	opt--structure that holds necessary input parameters

OUTPUTS:
	derivs--updated state vector

*********************************************************************************/

double * ebtel_derivs(double s[], double t, struct rk_params par, struct Option *opt)
{
 	//Declare variables
 	double p;
 	double n;
 	double T;
 	double rad;
 	double r3;
 	double f;
 	double f_eq;
 	double q;
	double pv;
 	double dpdt;
 	double dndt;
 	double dTdt;
	double *flux_ptr;
 	double *derivs = malloc(sizeof(double[3]));
 
 	//Unravel the state vector
 	p = s[0];
 	n = s[1];
 	T = s[2];
 	
 	//Compute the radiative loss function 
 	rad = ebtel_rad_loss(T,opt->rad_option);
 	
 	//Compute the coefficient r3
 	r3 = ebtel_calc_c1(t,T,n,par.L,rad,opt);
 	
 	//Compute heat flux
	flux_ptr = ebtel_calc_thermal_conduction(T,n,par.L,rad,r3,opt->sat_limit,opt->heat_flux_option);
	f = *(flux_ptr + 0);
	f_eq = *(flux_ptr + 1);
	free(flux_ptr);
	flux_ptr = NULL;
	
	//Set the heating depending on the time 
	q = ebtel_heating(t,opt);
	
	//Calculate the enthalpy flux
	pv = 0.4*(f_eq - f - par.flux_nt);
	
	//Now compute the derivatives of each of the quantities in our state vector
	dpdt = 	TWO_THIRDS*(q + (1. + 1./r3)*f_eq/par.L - (1. - 1.5*K_B*T/opt->energy_nt)*par.flux_nt/par.L);
	dndt = (pv*0.5/(par.r12*K_B*T*par.L) + par.flux_nt/opt->energy_nt/par.L);
	dTdt = T*(1./p*dpdt - 1./n*dndt);
	
	//Set the derivative state vector
	if(strcmp(opt->solver,"euler")==0)
	{
		derivs[0] = dpdt*(opt->tau) + p;
		derivs[1] = dndt*(opt->tau) + n;
		derivs[2] = derivs[0]/(2.0*derivs[1]*K_B);
	}
	else if(strcmp(opt->solver,"rka4") == 0 || strcmp(opt->solver,"rk4") == 0)
	{
		derivs[0] = dpdt;
		derivs[1] = dndt;
		derivs[2] = dTdt;
	}
	else
	{
		printf("Invalid solver option.\n");
		exit(0);
	}
	
	//Return the pointer
	return derivs;
}
 
  /**********************************************************************************
 
 FUNCTION NAME: ebtel_rk
 
 FUNCTION DESCRIPTION: This function implements a fourth-order Runge-Kutta routine in EBTEL. It will
 call the ebtel_rk_derivs function to compute derivatives of the necessary functions.
 This implementation is based on a routine written in MATLAB by A. L. Garcia (see 
 Numerical Methods for Physics, Prentice Hall, 1994).
  
 INPUTS:
	s--vector that stores the current state of the system
	n--number of variables in vector s
	t--current time
	tau--timestep
	par--structure that holds necessary parameters
	opt--structure that contains necessary input parameters	
 
 OUTPUTS:
 	s_out--updated state vector
	
 *********************************************************************************/
 
 double * ebtel_rk(double s[], int n, double t, double tau, struct rk_params par,struct Option *opt)
 {
 	//Declare variables
 	double half_tau;
 	double t_half;
 	double t_full;
 	double f_temp;
 	double *f1;
 	double *f2;
 	double *f3;
 	double *f4;
 	double s_temp[n];
 	double *s_out = malloc(sizeof(double[n]));
 	int i;
 	
 	//Set time variables
 	half_tau = 0.5*tau;
 	t_half = t + half_tau;
 	t_full = t + tau;
 	
 	//Compute the first function f1
 	f1 = ebtel_derivs(s,t,par,opt);
	//Make the temporary state vector
 	for(i=0; i<n; i++)
 	{
 		f_temp = *(f1 + i);
 		s_temp[i] = s[i] + half_tau*f_temp;
 	}
 	
 	//Compute the second function f2
 	f2 = ebtel_derivs(s_temp,t_half,par,opt);
	//Rebuild the temporary state vector
 	for(i=0; i<n; i++)
 	{
 		f_temp = *(f2 + i);
 		s_temp[i] = s[i] + half_tau*f_temp;
 	}
 	
 	//Compute the third function f3
 	f3 = ebtel_derivs(s_temp,t_half,par,opt);
	//Rebuild the temporary state vector
 	for(i=0; i<n; i++)
 	{
 		f_temp = *(f3 + i);
 		s_temp[i] = s[i] + half_tau*f_temp;
 	}
 	
 	//Compute the fourth function f4
 	f4 = ebtel_derivs(s_temp,t_full,par,opt);
	//Rebuild the temporary state vector
 	for(i=0; i<n; i++)
 	{
 		f_temp = *(f4 + i);
 		s_temp[i] = s[i] + tau*f_temp;
 	}
 	
 	//Now compute the final resulting state vector
 	for(i=0; i<n; i++)
 	{
 		s_out[i] = s[i] + 1./6.*tau*(*(f1+i) + *(f4+i) + 2.*( *(f2+i) + *(f3+i) ));
 	}
 	
 	//Free memory of all f functions
 	free(f1);
 	free(f2);
 	free(f3);
 	free(f4);
 	
 	//Return the final resulting state vector (pointer)
 	return s_out;
 	
 }
 
 /**********************************************************************************
 
 FUNCTION NAME: ebtel_rk_adapt
 
 FUNCTION DESCRIPTION: This function implements an adaptive step size for the Runge-Kutta
 routine used in EBTEL. It returns a structure pointer containing the current time, timestep,
 and state vector. It calls the ebtel_rk function. This implementation is based on a routine
 written in MATLAB by A. L. Garcia (see Numerical Methods for Physics, Prentice Hall, 1994). 
  
 INPUTS:
 	s--state vector
 	n--size of s
 	t--current time
 	tau--current time step
 	err--desired truncation error
 	par--structure used in computing derivatives
 	opt--structure containing input options
		
 OUTPUTS:
	rka_params--structure containing updated time, state, and timestep
	
*********************************************************************************/
 
 struct ebtel_rka_st *ebtel_rk_adapt(double s[], int n, double t, double tau, struct rk_params par, struct Option *opt)
 {
 	/**Declare variables**/
 	//Int
 	int i;
 	int j;
 	int max_try = 100;
	
 	//Double
 	double safe1 = 0.9;
 	double safe2 = 1.1;
	double safe3 = 4;
 	double scale;
 	double x_diff;
 	double error_ratio;
 	double t_save = t;
 	double time;
 	double half_tau;
 	double old_tau;
	double err = opt->rka_error;
 	double epsilon = 1.0e-16;
 	
	//Pointers
 	double *x_small_1;
 	double *x_small_2;
 	double *x_big;
 	
	//Arrays
 	double s_small_1[n];
 	double s_small_2[n];
 	double s_big[n];
 	
	//Structure
 	struct ebtel_rka_st *rka_params = malloc(sizeof(struct ebtel_rka_st));
 	//Reserve memory for structure members that will be set
 	rka_params->state=malloc(sizeof(double[n]));
 	
 	//Loop over maximum number of attempts to satisfy error bound
 	for(i=0; i<max_try; i++)
 	{
 		//Take two small steps
 		
 		//First small step
 		half_tau = 0.5*tau;
		
		x_small_1 = ebtel_rk(s,n,t_save,half_tau,par,opt);
		
 		//Unpack the x_small_1 pointer
 		for(j=0;j<n;j++)
 		{
 			s_small_1[j] = *(x_small_1 + j);
 		}
 		
 		//Update the time vector
 		time = t_save + half_tau;
 		
 		//Second small step
		x_small_2 = ebtel_rk(s_small_1,n,time,half_tau,par,opt);
		
		//Unpack the x_small_2 pointer
 		for(j=0;j<n;j++)
 		{
 			s_small_2[j] = *(x_small_2 + j);
 		}
 		
 		//Take single big step
 		x_big = ebtel_rk(s,n,t_save,tau,par,opt);
		
		//Unpack the x_big pointer
 		for(j=0;j<n;j++)
 		{
 			s_big[j] = *(x_big + j);
 		}
 		
 		//Update the time vector
 		time = t_save + tau;
 		
		error_ratio = 0.;
		//Compute estimated truncation error
		for(j=0;j<n;j++)
 		{
 			scale = err*(fabs(s_small_2[j]) + fabs(s_big[j]))/2.0;
			x_diff = s_small_2[j] - s_big[j];
 			//Return the maximum value of the error ratio
			error_ratio = ebtel_max_val(error_ratio,fabs(x_diff)/(scale + epsilon));
		}
		
 		//Estimate new tau value (including safety factors)
 		old_tau = tau;
 		tau = safe1*old_tau*pow(error_ratio,-1./5.);
 		tau = ebtel_max_val(tau,old_tau/safe2);	
		
 		//If error is acceptable, set our structure values and return the structure
 		if(error_ratio < 1)// && error_ratio != 0)
 		{
 			//Set the structure fields
			tau = ebtel_min_val(tau,safe3*old_tau);
			
 			rka_params->tau = tau;
 			for(j=0;j<n;j++)
 			{
 				rka_params->state[j] = s_small_2[j];
 			}
 			
 			//Free the memory of the pointers that were used
 			free(x_small_1);
 			x_small_1 = NULL;
 			free(x_small_2);
 			x_small_2 = NULL;
 			free(x_big);
 			x_big = NULL;
			
 			//Return the structure
 			return rka_params;
 		}
 		
 		//Free memory of the pointers that were used. They will be malloc'd on the next iteration
 		free(x_small_1);
 		x_small_1 = NULL;
 		free(x_small_2);
 		x_small_2 = NULL;
 		free(x_big);
 		x_big = NULL;
 	}
 	
 	//If we've finished the loop without meeting the error requirement, then return an error
 	printf("Error: Adaptive Runge-Kutta routine failed after %d iterations. Exiting the program\n",i);
	
	//Exit if the routine fails
	exit(0);
 	
 }
 