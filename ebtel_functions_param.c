/***********************************************************************************

FILENAME: ebtel_functions_param.c

AUTHOR Will Barnes

DATE: created: 7 March 2014

DESCRIPTION: This file contains all of the functions used to calculate important parameters
(c1,c2,c3,lambda) in EBTEL. All functions are described below. For additional details, 
see ebtel_main.c.

***********************************************************************************/

//Include appropriate header file
#include "ebtel_functions.h"

/***********************************************************************************

FUNCTION NAME: ebtel_calc_c1

FUNCTION_DESCRIPTION: This function calculates the ratio between radiative losses from the transition region and the corona The correspondence between the EBTEL code and Klimchuk et al.(2008) is c1=r3. Adapted from the function pro cacl_c1 in the original IDL EBTEL code. 

INPUTS:
	temp--temperature (K)
	den--electron number density (cm^-3).
	llength--loop half length (cm).
	rad--radiative loss function
	
OUTPUTS:
	c1--ratio of transition region to coronal radiative loss functions

***********************************************************************************/

double ebtel_calc_c1( double temp, double den, double llength, double rad )
{

	//Declare variables
	double sc;
	double f_eq_1;
	double r3_eqm_g;
	double r3_radn_g;
	double r3_eqm;
	double r3_radn;
	double n_eq_2;
	double noneq2;
	double r2;
	double r3;
	double r3_rad_0 = 0.6;	//radiative phase value, no gravity
	double r3_eqm_0 = 2.0;	//value in equilibrium with no gravity, -2/3 loss power law
	double l_fact_eq = 5.0;		//geometric factors for inclusion of gravitational effects
	double l_fact_rad = 5.0;	//l_fact^-1 = s/L*1/2 where <s/L> approx 0.4 so 1/l_fact apprrox 0.2 or 1/5
	
	//Calculate the scale height
	sc = ebtel_calc_lambda(temp);
	
	//Calculate r2 value
	r2 = ebtel_calc_c2();
	
	//Set equilibrium value of f
	f_eq_1 = -pow(den,2)*rad*llength;
	
	//Adjust values for gravity
	r3_eqm_g = r3_eqm_0*exp(4*sin(PI/l_fact_eq)*llength/(PI*sc));
	r3_radn_g = r3_rad_0*exp(4*sin(PI/l_fact_rad)*llength/(PI*sc));
	
	//Adjust for loss function
	int lossOff = 1;	//Use this to turn off (0) the loss function argument
	if (lossOff == 1)
	{
		r3_eqm = r3_eqm_g*1.95e-18/pow(temp,TWO_THIRDS)/rad;
		r3_radn = r3_radn_g*1.95e-18/pow(temp,TWO_THIRDS)/rad;
	}
	else
	{
		r3_eqm = r3_eqm_0;
		r3_radn = r3_rad_0;
	}
	
	//Calculate over/under density
	n_eq_2 = KAPPA_0*pow((temp/r2),SEVEN_HALVES)/(3.5*r3_eqm*rad*pow(llength,2));
	noneq2 = pow(den,2)/n_eq_2;
	
	//Use different values of r3 based on value of noneq2
	if (noneq2 < 1)
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

FUNCTION_DESCRIPTION: This function calculates the parameter c2 which is the ratio of average temperature to apex temperature. Adapted from the function pro calc_c2 in the original IDL EBTEL code.

INPUTS:
	
OUTPUTS:
	c2--ratio of tbar to ta

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

FUNCTION_DESCRIPTION: This function calculates the parameter c3 which is the ratio of base temperature to apex temperature. Adapted from the function pro calc_c3 in the original IDL EBTEL code.

INPUTS:
	
OUTPUTS:
	c3--ratio of t0 to ta

***********************************************************************************/

double ebtel_calc_c3(void)
{
	double c3;
	
	c3 = 0.6;
	//c3 = 0.5;	//Paper I value
	
	return c3;
	
}

/***********************************************************************************

FUNCTION NAME: ebtel_calc_lambda

FUNCTION_DESCRIPTION: This function calculates the scale height for a given temperature. Adapted from the function pro calc_lambda in the original IDL EBTEL code.

INPUTS:
	temp--temperature (K)
	
OUTPUTS:
	sc--scale height (cm)

***********************************************************************************/

double ebtel_calc_lambda( double temp )
{
	double sc;
	
	sc = (2*K_B*temp/M_P)/2.74e+4;
	
	return sc;
}