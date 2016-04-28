/***********************************************************************************

FILENAME: ebtel_functions_param.c

AUTHOR Will Barnes

DATE: created: 7 March 2014

***********************************************************************************/

//Include appropriate header file
#include "ebtel_functions.h"

/***********************************************************************************

FUNCTION NAME: ebtel_calc_c1

FUNCTION_DESCRIPTION: This function calculates the ratio between radiative losses 
from the transition region and the corona The correspondence between the EBTEL code 
and Klimchuk et al.(2008) is c1=r3. Adapted from the function pro cacl_c1 in the original
IDL EBTEL code. 

INPUTS:
	temp--temperature (K)
	den--electron number density (cm^-3).
	loop_length--loop half length (cm).
	rad--radiative loss function
	opt--structure that holds all inputs
	
OUTPUTS:
	c1--ratio of transition region to coronal radiative loss functions

***********************************************************************************/

double  ebtel_calc_c1( double t, double temp, double den, double loop_length, double rad, struct Option *opt )
{

	//Declare variables
	double sc;
	double r3_eqm_g;
	double r3_radn_g;
	double r3_eqm;
	double r3_radn;
	double n_eq_2;
	double noneq2;
	double r2;
	double r3;
	double r3_eqm_0;
	//double r3_rad_0 = 0.6;	//radiative phase value, no gravity
	//double r3_eqm_0 = 2.0;	//value in equilibrium with no gravity, -2/3 loss power law
	double l_fact_eq = 5.0;		//geometric factors for inclusion of gravitational effects
	double l_fact_rad = 5.0;	//l_fact^-1 = s/L*1/2 where <s/L> approx 0.4 so 1/l_fact apprrox 0.2 or 1/5
	
	//Calculate the scale height
	sc = ebtel_calc_lambda(temp);
	
	//Calculate r2 value
	r2 = ebtel_calc_c2();
	
	//Adjust values for sound speed correction
	if (strcmp(opt->r3_sound_speed_correction,"true")==0 || strcmp(opt->r3_sound_speed_correction,"True")==0)
	{
		double c_s,t_zero,tau_c1;
		int i;
		
		c_s = ebtel_calc_sound_speed(2.0*K_B*den*temp,den);
		tau_c1 = opt->tr_thickness*loop_length/c_s;
		
		t_zero = tau_c1;
		
		for(i=opt->num_events-1; i>=0; i--)
		{
			if(t>=opt->t_start_array[i])
			{
				t_zero = t - opt->t_start_array[i];
				break;
			}
		}
		
		if(t_zero < tau_c1)
		{
			r3_eqm_0 = (opt->r3_eqm_0a*(1.0-t_zero/tau_c1) + 2.0*opt->r3_eqm_0b*pow((t_zero/tau_c1),2))/(1.0 + pow((t_zero/tau_c1),2));
		}
		else
		{
			r3_eqm_0 = opt->r3_eqm_0b;
		}
	}
	else
	{
		r3_eqm_0 = opt->r3_eqm_0a;
	}
	
	//Adjust values for gravity
	if (strcmp(opt->r3_grav_correction,"true")==0 || strcmp(opt->r3_grav_correction,"True")==0)
	{
		r3_eqm_g = r3_eqm_0*exp(4*sin(PI/l_fact_eq)*loop_length/(PI*sc));
		r3_radn_g = opt->r3_rad_0*exp(4*sin(PI/l_fact_rad)*loop_length/(PI*sc));
	}
	else
	{
		r3_eqm_g = r3_eqm_0;
		r3_radn_g = opt->r3_rad_0;
	}
	
	//Adjust for loss function
	if (strcmp(opt->r3_loss_correction,"true")==0 || strcmp(opt->r3_loss_correction,"True")==0)
	{
		r3_eqm = r3_eqm_g*1.95e-18/pow(temp,TWO_THIRDS)/rad;
		r3_radn = r3_radn_g*1.95e-18/pow(temp,TWO_THIRDS)/rad;
	}
	else
	{
		r3_eqm = r3_eqm_g;
		r3_radn = r3_radn_g;
	}
	
	//Calculate over/under density
	n_eq_2 = KAPPA_0*pow((temp/r2),SEVEN_HALVES)/(SEVEN_HALVES*r3_eqm*rad*pow(loop_length,2));
	noneq2 = pow(den,2)/n_eq_2;
	
	//Use different values of r3 based on value of noneq2
	if (noneq2 < 1.0)
	{
		//Hot loops equilibrium value of c1
		r3 = r3_eqm;
	}
	else
	{	
		//Radiative loops transition from equilibrium (noneq2-->1)
		r3 = (2*r3_eqm + r3_radn*(noneq2-1))/(1+noneq2);
	}
	
	//Return the value of the parameter
	return r3;
	
}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_c2

FUNCTION_DESCRIPTION: This function calculates the parameter c2 which is the ratio of
average temperature to apex temperature. Adapted from the function pro calc_c2 in the
original IDL EBTEL code.

INPUTS:
	
OUTPUTS:
	c2--ratio of spatially averaged temperature to apex temperature

***********************************************************************************/

double ebtel_calc_c2( void )
{
	double c2;
	
	c2 = 0.9;
	//c2 = 0.87;	//Paper I value
	
	return c2;
}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_c3

FUNCTION_DESCRIPTION: This function calculates the parameter c3 which is the ratio of 
base temperature to apex temperature. Adapted from the function pro calc_c3 in the 
original IDL EBTEL code.

INPUTS:
	
OUTPUTS:
	c3--ratio of base temperature to apex temperature

***********************************************************************************/

double ebtel_calc_c3(void)
{
	double c3;
	
	c3 = 0.6;
	//c3 = 0.5;	//Paper I value
	
	return c3;
	
}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_sound_speed

FUNCTION_DESCRIPTION: This function calculates the sound speed in the coronal plasma.

INPUTS:
	p--pressure (dyne cm^-2)
	n--density (cm^-3)
	
OUTPUTS:
	cs--sound speed (cm s^-1)

***********************************************************************************/

double ebtel_calc_sound_speed(double p, double n)
{
	return pow(5.0/3.0*p/(M_P*n),0.5);
}


/***********************************************************************************

FUNCTION NAME: ebtel_calc_lambda

FUNCTION_DESCRIPTION: This function calculates the scale height for a given temperature.
Adapted from the function pro calc_lambda in the original IDL EBTEL code.

INPUTS:
	temp--temperature (K)
	
OUTPUTS:
	sc--scale height (cm)

***********************************************************************************/

double ebtel_calc_lambda( double temp )
{
	double sc;
	
	sc = (2.0*K_B*temp/M_P)/2.74e+4;
	
	return sc;
}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_abundance

FUNCTION_DESCRIPTION: This function calculates the effective Boltzmann constant and 
ion mass given some He/H abundance.

INPUTS:
	
OUTPUTS:

***********************************************************************************/

void ebtel_calc_abundance(void)
{
	double k_b = 1.38e-16;
	double m_p = 1.67e-24;
	
	//Calculate average ion mass
    double n_he_n_p = 0.075;   //He/p abundance.
    //double z_avg = (1.0 + 2.0*n_he_n_p)/(1.0 + n_he_n_p); //Include Helium
    double z_avg = 1.; //For Hydrad comparison.
    double kb_fact = 0.5*(1.0+1.0/z_avg);
    //double m_fact = (1.0 + n_he_n_p*4.0)/(2.0 + 3.0*n_he_n_p); //Include Helium
    double m_fact = (1 + n_he_n_p*4.)/2.; //For Hydrad comparison
	
	//Set global variables
    M_P = m_p*m_fact*(1.0 + z_avg)/z_avg; //Average ion mass
	K_B = k_b*kb_fact; //Modify equation of state for non-e-p plasma
}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_ic

FUNCTION_DESCRIPTION: This function calculates the initial temperature using the static equilibrium
conditions for the hydrodynamic equations. 

INPUTS:
	kpar--temperature bins used to calculate the radiative loss function
	r3--ratio of TR and coronal radiative losses (coefficient c1 from Papers I,II)
	loop_length--loop half-length (Mm)
	opt--Option structure of input values
	
OUTPUTS:
	return_array--array that holds r3, rad, t, n, p, and v

***********************************************************************************/

double * ebtel_calc_ic(double r3, double loop_length, struct Option *opt)
{
	//Variable declarations for both cases
	double *return_array = malloc(sizeof(double[6]));
	double r2 = ebtel_calc_c2();
	double heat = ebtel_heating(0,opt);
	
	if(strcmp(opt->ic_mode,"force") == 0 || strcmp(opt->ic_mode,"st_eq") == 0)
	{
		//Variable declarations
		int i;
		double tt_old;
		double nn_old;
		double tt_new;
		double rad;
		double nn;
		double p;
		double v;
		double err;
		double err_n;
		double tol;
	
		//Check if the heating array begins with a zero. If so, return an error.
		if (heat == 0.)
		{
			printf("ERROR! No initial loop heating: heat(0)=0. Provide valid heating input. Exiting the program\n");
			exit(0);
		}

		//First set up trial values for static equilibrium (i.e. d/dt = 0)
		tt_old = r2*pow(3.5*r3/(1 + r3)*loop_length*loop_length*heat/KAPPA_0,TWO_SEVENTHS);
		rad = ebtel_rad_loss(tt_old,opt->rad_option);
		nn = pow(heat/((1+r3)*rad),0.5);
		nn_old = nn;

		//Compute initial values for parameters t and n by iterating on temperature (tt) and R_tr/R_c (r3)
		tol = 1e+3;		//error tolerance

		for(i=0; i<=100; i++)
		{
			r3 = ebtel_calc_c1(0.0,tt_old,nn,loop_length,rad, opt);										//recalculate r3 coefficient
			tt_new = r2*pow((3.5*r3/(1+r3)*pow(loop_length,2)*heat/KAPPA_0),TWO_SEVENTHS);		//temperature at new r3
			rad = ebtel_rad_loss(tt_new,opt->rad_option);									//radiative loss at new temperature
			nn = pow(heat/((1+r3)*rad),0.5);													//density at new r3 and new rad
			err = tt_new - tt_old;																//difference between t_i, T_i-1
			err_n = nn - nn_old;	
			//Break the loop if the error gets below a certain threshold
			if(fabs(err)<tol && fabs(err_n)<tol)
			{
				//Set parameters and break loop
				tt_old = tt_new;
				nn_old = nn;
				break;
			}
			tt_old = tt_new;
			nn_old = nn;
		}
	
		//Calculate the density
		nn = pow(heat/((1+r3)*rad),0.5);
		
		//To use parameters from configuration file
		if(strcmp(opt->ic_mode,"force") == 0)
		{
			tt_old = opt->T0;
			nn = opt->n0;
		}
		
		//Calculate resulting pressure, velocity
		p = 2*K_B*nn*tt_old;
		v = 0;
	
		//Set array values
		return_array[0] = r3;
		return_array[1] = rad;
		return_array[2] = tt_old;
		return_array[3] = nn;
		return_array[4] = p;
		return_array[5] = v;
	}
	else if(strcmp(opt->ic_mode,"scaling") == 0)
	{
		//Variable declarations
		double lambda_0;
		double bb;
		double t_0;
		double p_0;
		double n_0;
		double v_0;
		double rad;
		
		//Alternatively, we could use the scaling laws to determine our initial conditions
		lambda_0 = 1.95e-18;			//lambda = lambda_0*T
		bb = -TWO_THIRDS;//-0.5			//power law for radiative loss function
		t_0 = r2*pow((3.5/KAPPA_0*heat),TWO_SEVENTHS)*pow(loop_length,2.0*TWO_SEVENTHS);
		p_0 = pow(r2,-SEVEN_HALVES*0.5)*pow(8.0/7.0*KAPPA_0/lambda_0,0.5)*K_B*pow(t_0,((11.0-2.0*bb)/4.0))/loop_length;
		n_0 = 0.5*p_0/(K_B*t_0);
		v_0 = 0;
		
		//Set array values
		rad = ebtel_rad_loss(t_0,opt->rad_option);
		return_array[0] = ebtel_calc_c1(0.0,t_0,n_0,loop_length,rad,opt);
		return_array[1] = rad;
		return_array[2] = t_0;
		return_array[3] = n_0;
		return_array[4] = p_0;
		return_array[5] = v_0;
	}
	
	
	return return_array;
}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_thermal_conduction

FUNCTION_DESCRIPTION: This function calculates the thermal conduction, commonly denoted
F_c. It uses the classical Spitzer-Harm approximation and, when specified, uses a flux limit
approximation to prevent run-away cooling in situations where the density is not sufficient
to support such a flux.

INPUTS:
	T--temperature (K)
	n--density (cm^-3)
	loop_length--loop half-length (cm)
	rad--radiative loss
	r3--the c1 ebtel parameter (ratio of TR to coronal losses)
	heat_flux_option--choose to use 'classical' or 'dynamic' heat flux calculation
	
OUTPUTS:
	flux_ptr--pointer to array that holds heat flux and equilibrium heat flux

***********************************************************************************/

double * ebtel_calc_thermal_conduction(double T, double n, double loop_length, double rad, double r3, double sat_limit, char *heat_flux_option)
{
	
	//Declare variables
	double c1;
	double c_sat;
	double f_cl;
	double f_sat;
	double f;
	double f_eq;
	double r2 = ebtel_calc_c2();
	double *flux_ptr = malloc(sizeof(double[2]));
	
	//Set up thermal conduction parameters (NEED TO CHANGE FOR e- AND ion)
	c1 = -TWO_SEVENTHS*KAPPA_0;
	c_sat = -1.5*pow(K_B,1.5)/pow(M_EL,0.5);
	
	//Set up thermal conduction at the base
	f_cl = c1*pow(T/r2,SEVEN_HALVES)/loop_length;	//Classical heat flux calculation
	
	//Decide on whether to use classical or dynamic heat flux
	if(strcmp(heat_flux_option,"classical")==0)
	{
		f = f_cl;
	}
	else if(strcmp(heat_flux_option,"limited")==0)
	{
		//Compute flux limit
		f_sat = sat_limit*c_sat*n*pow(T,1.5);
		
		//Compute final flux value		
		f= -f_cl*f_sat/pow((pow(f_cl,2.) + pow(f_sat,2)),0.5);		
	}
	else
	{
		printf("Invalid heat flux option.\n");
		exit(0);
	}
	
	//Calculate equilibrium thermal conduction at base (-R_tr in Paper I)
	f_eq = -r3*pow(n,2.)*rad*loop_length;
	
	//Set the flux array
	flux_ptr[0] = f;
	flux_ptr[1] = f_eq;
	
	//Return the array
	return flux_ptr;
	
}

