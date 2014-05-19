/***********************************************************************************

FILENAME: ebtel_functions_loop.c

AUTHOR Will Barnes

DATE: created: 7 March 2014

DESCRIPTION: This file contains all of the major functions used in the EBTEL model.
Additional functions can be found in ebtel_functions_util.c and ebtel_functions_param.c.
Full function descriptions can be found below. See ebtel_main.c for additional details.
***********************************************************************************/

//Include appropriate header file
#include "ebtel_functions.h"

/***********************************************************************************

FUNCTION NAME: ebtel_loop_solver

FUNCTION DESCRIPTION: Given some heating input, this function solves the 0-D hydrodynamic equations spatially averaged over a loop half-length, computing the temperature, density, and pressure as functions of time. Given appropriate initial conditions (see below), the differential emission measure of the transition region and corona can also be calculated. 

INPUTS:
   ttime--time array
   heat--heating rate array (erg cm^-3 s^-1)
   looplength--loop half length (cm) (top of the chromosphere to apex)
   opt--keyword structure: (set to ON (1) or OFF (0))
        classical--set to use the unsaturated classical heat flux
        dynamic--set to use dynamical r1 and r2
        dem_old--set to use old technique of computing DEM(T) in the TR
        rtv--set to use Rosner, Tucker, & Vaiana radiative loss function
        usage--this keyword is not optional and must be set to determine
               proper usage. See USAGE section below.
        Following value is not a keyword and contains actual values:
        energy_nt--mean energy of the nonthermal electrons (keV)

OUTPUTS:
   t--temperature (K) array corresponding to time array (avg. over coronal
      section of the loop)
   n--electron number density array (cm^-3) (coronal average)
   p--pressure array (dyn cm^-2) (coronal avg.)
   v--velocity array (cm s^-1) (r4 * velocity at base of corona)
   dem_tr--differential emission measure of transition region, dem(time,T), both legs
             (dem = n^2 * ds/dT  cm^-5 K^-1)
             (Note:  dem_tr is not reliable when a nonthermal electron flux is used.)
   dem_cor--differential emission measure of corona, dem(time,T), both legs
             (Int{dem_cor+dem_tr dT} gives the total emission measure of a 
             loop strand having a cross sectional area of 1 cm^2)
   logtdem--logT array corresponding to dem_tr and dem_cor
   f_ratio--ratio of heat flux to equilibrium heat flux
             (ratio of heat flux to tr. reg. radiative loss rate)
   rad_ratio--ratio of tr. reg. radiative loss rate from dem_tr and from r3*(coronal rate)
   cond--conductive loss from the corona
   rad_cor--coronal radiative loss
***********************************************************************************/

struct ebtel_params_st *ebtel_loop_solver( int ntot, double loop_length, double total_time, double time[], double heat[], struct Option opt)
{
	/***********************************************************************************
								Variable Declarations
	***********************************************************************************/
	
	/*****Variable declarations*****/
	//Int
	int nk;
	int i;	//index over ntot
	int j;	//index over 451
	int flag_dem_tr;
	int j_min;
	int j_max;
	
	//Pointers
	double *kptr;
	double *state_ptr;
	double *log_tdem_ptr;
	
	//Double
	double r1;
	double r1_tr;
	double r2;
	double r3;
	double r4;
	double rad;
	double c1;
	double c_sat;
	double root_c2;
	double c3;
	double c4;
	double tt_old;
	double tt_new;
	double tt;
	double nn;
	double nn_old;
	double err;
	double err_n;
	double tol;
	double sc;
	double lambda_0;
	double bb;
	double q_0;
	double t_0;
	double v_0;
	double p_0;
	double n_0;
	double f_cl;
	double f;
	double f_sat;
	double sat_limit;
	double r12;
	double r12_tr;
	double f_eq;
	double pv;
	double cf;
	double t_max;
	double t_min;
	double em;
	double dem0;
	double delta_t;
	double rad_loss;
	double t;
	double n;
	double ta;
	double na;
	double p;
	double pa;
	double v;
	double cond;
	double rad_cor;
	double rad_ratio;
	double f_ratio;
	double tau;
	
	//Double array (single index)
	double c_array[3];
	double f_array[3];
	double state[3];
	double log_tdem[451];
	double tdem[451];
	double root_tdem[451];
	double fourth_tdem[451];
	double rad_dem[451];
	double root_rad_dem[451];
	double dem_cor_minus[ntot];
	double dem_tr_minus[ntot];
	double dem_tot_minus[ntot];
	
	//Double array (dynamically allocate memory to avoid segmentation fault)
	double **dem_tr;
	dem_tr = malloc(ntot*sizeof(double *));
	double **dem_cor;
	dem_cor = malloc(ntot*sizeof(double *));
	for(i = 0; i<ntot; i++)
	{
		dem_tr[i] = malloc(451*sizeof(double));
		dem_cor[i] = malloc(451*sizeof(double));
	}
	
	//struct
	struct rk_params par;
	struct ebtel_params_st *param_setter = malloc(sizeof(struct ebtel_params_st));
	assert(param_setter != NULL);
	
	//Reserve memory for structure members that will be set
	param_setter->time = malloc(sizeof(double[ntot]));
	param_setter->heat = malloc(sizeof(double[ntot]));
	param_setter->temp = malloc(sizeof(double[ntot]));
	param_setter->ndens = malloc(sizeof(double[ntot]));
	param_setter->press = malloc(sizeof(double[ntot]));
	param_setter->vel = malloc(sizeof(double[ntot]));
	param_setter->tapex = malloc(sizeof(double[ntot]));
	param_setter->napex = malloc(sizeof(double[ntot]));
	param_setter->papex = malloc(sizeof(double[ntot]));
	param_setter->coeff_1 = malloc(sizeof(double[ntot]));
	param_setter->cond = malloc(sizeof(double[ntot]));
	param_setter->rad_cor = malloc(sizeof(double[ntot]));
	
	if(opt.usage == 4 || opt.usage == 1)
	{
		if(opt.usage == 4)
		{	
			param_setter->f_ratio = malloc(sizeof(double[ntot]));
			param_setter->rad_ratio = malloc(sizeof(double[ntot]));
		}
		param_setter->logtdem = malloc(sizeof(double[451]));
		param_setter->dem_tr_log10mean = malloc(sizeof(double[451]));
		param_setter->dem_cor_log10mean = malloc(sizeof(double[451]));
		param_setter->dem_tot_log10mean = malloc(sizeof(double[451]));
	}
	
	/***********************************************************************************
									Initial Parameters
	***********************************************************************************/	
	
	//Calculate initial r values. In Klimchuk et al. (2008), these are denoted by c1,c2,c3. 
	//See ebtel_main.c documentation for correspondence between variables.
	r1 = ebtel_calc_c3();	//ratio of base to apex temperature; c3 in Klimchuk et al. (2008), Paper I from here on.
	r1_tr = r1;				//ratio of temperature at top of TR to apex temperature
	r2 = ebtel_calc_c2();	//ratio of average to apex temperature; c2 in Paper I
	r3 = 2;					//ratio of TR to coronal radiative losses; c1 in Paper I
	r4 = 1.0;				//ratio of average to base velocity
	
	//Set temperature bins.
	//Lengths of the kpar array are different depending on the loss function we use.
	if (opt.rtv== 0)
	{nk = 7;
	}
	else
	{nk = 6;
	}
	//Declare kpar array which will set the bins based on our loss function choice
	double kpar[nk];
	//Call the ebtel_kpar_set function and pass the pointer kptr
	kptr = ebtel_kpar_set(opt.rtv);
	//Set the kpar array;
	for(i = 0; i < nk; i++)
	{
		kpar[i] = *(kptr + i);
	}
	//Calculate radiative loss function.
	rad = ebtel_rad_loss(1e+6,kpar,opt.rtv);
	
	//Set up thermal conduction parameters
	c1 = -TWO_SEVENTHS*KAPPA_0;
	c_sat = -1.5*pow(K_B,1.5)/pow(9.1e-28,0.5);
	//sat_limit = 0.1667;
	sat_limit = 1;	//HYDRAD value
	
	/***********************************************************************************
						Set up DEM in Transition Region
	***********************************************************************************/
	if (opt.usage == 1 || opt.usage == 4)
	{	
		//Define temperature arrays for plotting and calculating DEM
		log_tdem_ptr = ebtel_linspace(0,450,451);		//set the linspace pointer using the ebtel_linspace function
		for(i = 0; i<451; i++)
		{
			log_tdem[i] = *(log_tdem_ptr + i)/100 + 4;	//log of T in the region in which DEM is defined
			param_setter->logtdem[i] = log_tdem[i];		//Save logtdem to our parameter structure
			tdem[i] = pow(10,log_tdem[i]);				//T in the region in which DEM is defined
			root_tdem[i] = sqrt(tdem[i]);				//T^1/2 in the region in which DEM is defined
			fourth_tdem[i] =  pow(tdem[i],0.25);			//T^1/4 in the region in which DEM is defined
		}
		
		//These coefficients will be used in the old method of calculating DEM
		root_c2 = pow((KAPPA_0/(20*K_B)),0.5)/K_B;		//Calculate the root of c2 to avoid overflow in our calculation of dem_ev
		c3 = -5*K_B;
		c4 = pow((KAPPA_0/14),0.5)/K_B;
		//Make c_array for use in ebtel_calc_tr_dem function
		c_array[0] = root_c2; c_array[1] = c3; c_array[2] = c4;
		
		//Radiation in the transition region. This loop just calculates the radiative loss function in the TR
		for(i = 0; i<451; i++)
		{
			rad = ebtel_rad_loss(tdem[i],kpar,opt.rtv);		//Calculate the radiative loss function for temperature tdem[i]
			rad_dem[i] = rad;								//Set radiative loss function in the TR
			if (tdem[i] < 1e+4)
			{
				rad_dem[i] = 1;								//Check to see if we are outside the allowed temperature range
			}
			root_rad_dem[i] = sqrt(rad_dem[i]);
		}
	}
	
	/***********************************************************************************
							Initial Static Equilibrium
	***********************************************************************************/
	/*
	2 possible methods: (A) use EBTEL equilibrium (recommended) (B) use scaling laws
	*/
	
	//Check if the heating array begins with a zero. If so, return an error.
	if (heat[0] == 0)
	{
		printf("No initial loop heating: heat(0)=0. Provide valid heating input.\n");
	}

	//(A) EBTEL Equilibrium

	//First set up trial values for static equilibrium (i.e. d/dt = 0)
	tt_old = r2*pow(3.5*r3/(1 + r3)*loop_length*loop_length*heat[0]/KAPPA_0,TWO_SEVENTHS);
	printf("tt_old = %e\n",tt_old);
	rad = ebtel_rad_loss(tt_old,kpar,opt.rtv);
	nn = pow(heat[0]/((1+r3)*rad),0.5);
	nn_old = nn;
	
	/*
	//DEBUG
	printf("Print seed values to debug iteration on tt\n");
	printf("tt_old = %e\n",tt_old);
	printf("rad = %e\n",rad);
	printf("nn_old = %e\n",nn_old);
	*/

	//Compute initial values for parameters t and n by iterating on temperature (tt) and R_tr/R_c (r3)
	tol = 1e+3;		//error tolerance

	for(i=0; i<=100; i++)
	{
		r3 = ebtel_calc_c1(tt_old,nn,loop_length,rad);										//recalculate r3 coefficient
		tt_new = r2*pow((3.5*r3/(1+r3)*pow(loop_length,2)*heat[0]/KAPPA_0),TWO_SEVENTHS);	//temperature at new r3
		rad = ebtel_rad_loss(tt_new,kpar,opt.rtv);											//radiative loss at new temperature
		nn = pow(heat[0]/((1+r3)*rad),0.5);												//density at new r3 and new rad
		err = tt_new - tt_old;															//difference between t_i, T_i-1
		err_n = nn - nn_old;	
		//Break the loop if the error gets below a certain threshold
		if(fabs(err)<tol)// && fabs(err_n)<tol)
		{
			printf("r3 = %e\n",r3);													//display calculated parameters
			printf("tt_new = %e\n",tt_new);
			printf("tt_old = %e\n",tt_old);
			printf("err = %e\n",err);
			printf("Broke on iteration %d\n",i);
			tt_old = tt_new;
			nn_old = nn;
			break;
		}
		tt_old = tt_new;
		nn_old = nn;
	}

	//To use parameters consistent with the cases invoked in Paper II, we read in initial values for n,T rather than
	//calculating them using scaling laws or static equilibrium
	if(opt.mode==1)
	{
		tt = opt.T0;
		nn = opt.n0;
	}
	else
	{
		tt = tt_old;
		nn = pow(heat[0]/((1+r3)*rad),0.5);
	}
	
	//Print out the coefficients that we are starting the model with
	printf("********************************************************************\n");
	printf("Model Parameters\n");
	printf("r1 = %e\n",r1);
	printf("r2 = %e\n",r2);
	printf("r3 = %e\n",r3);
	printf("********************************************************************\n");
	
	//Set some initial parameters before iterating in time
	t = tt;
	n = nn;
	p = 2*K_B*n*t;
	v = 0;
	ta = t/r2;
	sc = ebtel_calc_lambda(t);
	na = n*r2*exp(-2.0*loop_length/(PI*sc)*(1.0-sin(PI/5.0)));
	pa = 2*K_B*na*ta;
	
	//Set the initial values of our parameter structure
	param_setter->temp[0] = t;
	param_setter->ndens[0] = n;
	param_setter->press[0] = p;
	param_setter->vel[0] = v;
	param_setter->tapex[0] = ta;
	param_setter->napex[0] = na;
	param_setter->papex[0] = pa;
	param_setter->coeff_1[0] = r3;
	param_setter->heat[0] = heat[0];
	param_setter->time[0] = 0;
	
	//Print out the coefficients that we are starting the model with
	printf("********************************************************************\n");
	printf("Model Parameters\n");
	printf("L = %e\n",loop_length);
	printf("Q = %e\n",heat[0]);
	printf("T, Ta = %e, %e\n",param_setter->temp[0],param_setter->tapex[0]);
	printf("n, na = %e, %e\n",param_setter->ndens[0],param_setter->napex[0]);
	printf("********************************************************************\n");
	
	//Alternatively, we could use the scaling laws to determine our initial conditions
	lambda_0 = 1.95e-18;			//lambda = lambda_0*T
	bb = -TWO_THIRDS;//-0.5				//power law for radiative loss function
	q_0 = heat[0];
	t_0 = r2*pow((3.5/KAPPA_0*heat[0]),TWO_SEVENTHS)*pow(loop_length,2.0*TWO_SEVENTHS);
	p_0 = pow(r2,-SEVEN_HALVES*0.5)*pow(8.0/7.0*KAPPA_0/lambda_0,0.5)*K_B*pow(t_0,((11.0-2.0*bb)/4.0))/loop_length;
	n_0 = 0.5*p_0/(K_B*t_0);
	v_0 = 0;
	
	//Print scaling law values to the screen
	printf("********************************************************************\n");
	printf("Scaling Law Values\n");
	printf("T_0 = %e\n",t_0);
	printf("P_0 = %e\n",p_0);
	printf("n_0 = %e\n",n_0);
	printf("********************************************************************\n");
	
	/***********************************************************************************
							Time-dependent Heating
	***********************************************************************************/
	
	//Set the structure members of the par structure. This is used in both the RK and Euler methods.
	par.L = loop_length;
	par.kpar = kptr;
	par.r12 = r1/r2;
	par.r2 = r2;
	par.sat_limit = sat_limit;
	par.c_sat = c_sat;
	par.c1 = c1;
	
	//Set the initial timestep from opt structure.
	//This will be static if we are not using our adaptive solver
	tau = opt.tau;
	
	//Begin the loop over the timesteps
	for(i = 0; i < ntot-1; i++)
	{
		//Update the parameter structure
		par.q1 = heat[i+1];
		par.q2 = heat[i];
		
		//Save time and heat to main data structure
		param_setter->heat[i+1] = heat[i+1];
		param_setter->time[i+1] = time[i+1];
		
		//Set up non-thermal electron flux for usage option 3
		if(opt.usage==3)
		{
			par.flux_nt = loop_length*heat[i]/10;
		}
		else
		{
			par.flux_nt = 0;
		}
		
		//Set up thermal conduction at the base
		f_cl = c1*pow(t/r2,SEVEN_HALVES)/loop_length;	//Classical heat flux calculation

		//Decide on whether to use classical or dynamic heat flux
		if(opt.classical==1)
		{
			f = f_cl;
		}
		else
		{
			f_sat = sat_limit*c_sat*n*pow(t,1.5);
			f = -f_cl*f_sat/pow((pow(f_cl,2) + pow(f_sat,2)),0.5);
		}
		par.f = f;

		//Calculate radiative loss
		rad = ebtel_rad_loss(t,kpar,opt.rtv);

		//Calculate coefficients r1, r2, r3 (c3, c2, c1)
		r3 = ebtel_calc_c1(t,n,loop_length,rad);
		par.r3 = r3;
		param_setter->coeff_1[i+1] = r3;
		r2 = ebtel_calc_c2();
		r1 = ebtel_calc_c3();
		r12 = r1/r2;
		r12_tr = r1_tr/r2;

		//Calculate equilibrium thermal conduction at base (-R_tr in Paper I)
		f_eq = -r3*n*n*rad*loop_length;
		par.f_eq = f_eq;
		
		//Calculate pv quantity to be used in velocity calculation
		pv = 0.4*(f_eq - f - par.flux_nt);
		
		/*****Step parameters forward in time using (1) RK method or (0) Euler stepper method*****/
		
		//Update the state vector
		state[0] = p;
		state[1] = n;
		state[2] = t;
		
		if(opt.solver==0)	//Euler solver
		{	
			//Call the Euler routine
			state_ptr = ebtel_euler(state,tau,par,opt);
		}
		else //if(opt.solver==1)	//RK routine
		{	
			//Call the RK routine
			state_ptr = ebtel_rk(state,3,time[i],tau,par,opt);	
		}

		//Update p,n,t and save to structure
		p = *(state_ptr + 0);
		param_setter->press[i+1] = p;
		n = *(state_ptr + 1);
		param_setter->ndens[i+1] = n;
		t = *(state_ptr + 2);
		param_setter->temp[i+1] = t;
		
		//Free the state ptr as it will be malloc'd again on the next go around
		free(state_ptr);
		state_ptr = NULL;
		
		v = pv/p; 			//calculate new velocity
		param_setter->vel[i+1] = v*r4;
		
		//Calculate new scale height
		sc = ebtel_calc_lambda(t);
		
		//Calculate apex quantities
		ta = t/r2;
		param_setter->tapex[i+1] = ta;
		na = n*r2*exp(-2.0*loop_length*(1.0-sin(PI/5.0))/(PI*sc));
		param_setter->napex[i+1] = na;
		pa = 2*K_B*na*ta;
		param_setter->papex[i+1] = pa;
		
		/*****Differential Emission Measure Calculation*****/
		//Check usage variable to determine whether we are calculating TR DEM
		if(opt.usage == 1 || opt.usage == 4)
		{	
			//Transition region
			if(r12_tr*t > tdem[450])
			{
				printf(" Transition region T = %e K outside of DEM range\n",r12_tr*t);
			}
			
			if(f != f_eq)
			{
				cf = f*f_eq/(f - f_eq);
			}
			else
			{
				cf = 1e+10*f;
			}
			
			//Make f array for ebtel_calc_tr_dem function
			f_array[0] = f;
			f_array[1] = f_eq;
			f_array[2] = cf;
			
			//Initialize dem_tr flag to zero. We want to know if dem_tr takes on negative values and if so we need to reset them.
			flag_dem_tr = 0;
			
			for(j=0; j<451; j++)
			{
				//Check to see whether we are in the TR. If so, calculate dem_TR. Note: r12_tr*t[i] = T_0
				if( tdem[j] < r12_tr*t )
				{
					//Make call to function that calculates DEM for TR using method specified by opt.dem_old.
					dem_tr[i][j] = ebtel_calc_tr_dem(tdem[j],n,v,p,loop_length,sc,rad_dem[j],c_array,f_array,opt.dem_old);
					
					//Check whether the dem is less than zero and set the flag if the condition holds
					if(dem_tr[i][j] < 0)
					{
						flag_dem_tr = 1;
						printf("***** Negative DEM at i = %d\n", i);
						break;
					}
					
					//Check whether it is classical or dynamic to set the f_ratio value
					if(opt.classical == 1)
					{
						f_ratio = f/f_eq;
					}
					else
					{
						f_ratio = f_cl/f_sat;
					}
				}
			}
			
			//If we tripped the flag_dem_tr, then we need to reset dem_tr
			if(flag_dem_tr == 1)
			{
				for(j=0; j<451; j++)
				{
					if(i != 0)
					{
						dem_tr[i][j] = dem_tr[i-1][j];
					}
					else
					{
						dem_tr[i][j] = 0;
					}
				}
			}
			
			//Corona (EM distributed uniformly over temperature interval [tmin,tmax])
			t_max = ebtel_max_val(t/r2,1e+4);
			t_min = ebtel_max_val(t*(2.0 - 1/r2),1e+4);
			j_max = (log10(t_max) - 4.0)*100;
			j_min = (log10(t_min) - 4.0)*100;
			
			//Calculate total emission measure
			em = 2*pow(n,2)*loop_length;			//factor of 2 for both legs
			
			delta_t = 1e4*(pow(10.0,((j_max + 0.5)/100.0)) - pow(10.0,((j_min - 0.5)/100.0)));
        	dem0 = em/delta_t;
        	
        	//Set every value between j_min and j_max DEM in the corona to dem0
        	for(j=j_min; j<=j_max; j++)
        	{
        		dem_cor[i][j] = dem0;
        	}
        	
        	//Transition region radiation losses based on DEM
        	if(opt.usage == 4)
        	{
        		rad_loss = 0;		//Initialize rad_loss to 0
        		
        		//Sum up radiative losses in the transition region
        		for(j=0; j<451; j++)
        		{
        			if(tdem[j] < r12_tr*t)
        			{
        				rad_loss += dem_tr[i][j]*rad_dem[j]*tdem[j]*0.01*log(10.0);
        			}
        		}
        		
        		//Compute rad_ratio
        		rad_ratio = -rad_loss/f_eq;
        		param_setter->rad_ratio[i] = rad_ratio;
        		param_setter->f_ratio[i] = f_ratio;
        	}
		}
		
		//Set the conductive loss from the corona 
		cond = f;
		param_setter->cond[i] = f;
		//Set the coronal radiative loss value
		rad_cor = f_eq/r3;
		param_setter->rad_cor[i] = rad_cor;
		
		//Increment the counter
		//i++;
	}
	
	//End of loop
	
	//Need to set final values for transition region and coronal radiative losses
	if(opt.usage == 1 || opt.usage == 4)
	{
		for(j = 0; j < 451; j++)
		{
			dem_tr[ntot-1][j] = dem_tr[ntot-2][j];
			dem_cor[ntot-1][j] = dem_cor[ntot-2][j];
			
			//Create a single dimensional array from a doubly indexed array
			for(i = 0; i<ntot; i++)
			{
				dem_cor_minus[i] = dem_cor[i][j];
				dem_tr_minus[i] = dem_tr[i][j];
				dem_tot_minus[i] = dem_tr[i][j] + dem_cor[i][j];
			}
			
			//Compute the mean of each of our newly created arrays over the index i, compute the log10, and store it as an entry in the array dem_log10mean_{region}. We will return these arrays 
			param_setter->dem_cor_log10mean[j] = log10(ebtel_avg_val(dem_cor_minus,ntot));
			param_setter->dem_tr_log10mean[j] = log10(ebtel_avg_val(dem_tr_minus,ntot));
			param_setter->dem_tot_log10mean[j] = log10(ebtel_avg_val(dem_tot_minus,ntot));
		}
	}
	
	//Free up memory used by ebtel_kpar_set and ebtel_linspace functions
	free(kptr);
	kptr = NULL;
	for(i=0;i<ntot;i++)
	{
		free(dem_tr[i]);
		dem_tr[i] = NULL;
		free(dem_cor[i]);
		dem_cor[i] = NULL;
	}
	free(dem_tr);
	dem_tr = NULL;
	free(dem_cor);
	dem_cor = NULL;
	if(opt.usage==1 || opt.usage==4)
	{
		free(log_tdem_ptr);
		log_tdem_ptr = NULL;
	}
	
	//Exit and return the structure that has been set appropriately
	return param_setter;
}

/***********************************************************************************

FUNCTION NAME: ebtel_kpar_set

FUNCTION DESCRIPTION: This function sets the kpar array to be used later in calculations of the radiative loss function. These are essentially the temperature bins used in calculating the radiative loss function. Either the RK or RTV method is used.

INPUTS:
	rtv_opt--use either the RK or RTV method 
OUTPUTS:
	*kpar--pointer for the kpar array

***********************************************************************************/

double * ebtel_kpar_set(int rtv_opt)
{	
	//Check option input to decide which method to use
	if (rtv_opt == 0)
	{	
		double *kpar = malloc(sizeof(double[7]));
		//Raymond-Klimchuk Loss function
		kpar[0] = 1.e+4;
		kpar[1] = 9.3325e4;
		kpar[2] = 4.67735e5;
        kpar[3] = 1.51356e6;
        kpar[4] = 3.54813e6;
        kpar[5] = 7.94328e6;
        kpar[6] = 4.28048e7;
        
        return kpar;
	}
	else
	{
		double *kpar = malloc(sizeof(double[6]));
		//Rosner-Tucker-Vaiana Loss function
		kpar[0] = pow(10,4.3);
        kpar[1] = pow(10,4.6);
        kpar[2] = pow(10,4.9);
        kpar[3] = pow(10,5.4);
        kpar[4] = pow(10,5.75);
        kpar[5] = pow(10,6.3);
        
        return kpar;
	}
}

/***********************************************************************************

FUNCTION NAME: ebtel_rad_loss

FUNCTION_DESCRIPTION: This function calculates the radiative loss function using either the Raymond-Klimchuk (RK) method or the Rosner-Tucker-Vaiana (RTV) method. Adopted from pro radLoss from the original EBTEL IDL code.

INPUTS:
	temp--temperature (K)
	kpar--holds the kpar structure if it has been previously set; placeholder otherwise.
	rtv_opt--option to use the RTV method (1) or the RK method (1).
	
OUTPUTS:
	rad--radiative loss function

***********************************************************************************/

double ebtel_rad_loss( double temp, double kpar[], int rtv_opt)
{
	//Declare rad 
	double rad;
	
	//Find out which method is being used
	if (rtv_opt == 0)
	{
		//RK loss function
    	if ( temp > kpar[6] ){ 
        	rad = 1.96e-27*sqrt(temp);
        }
    	else if ( temp > kpar[5] ){ 
        	rad = 5.4883e-16/temp;
        }
    	else if ( temp > kpar[4] ){
        	rad = 3.4629e-25*(pow(temp,0.333));
        }
    	else if ( temp > kpar[3] ){ 
        	rad = 3.5300e-13/(pow(temp,1.5));
        }
    	else if ( temp > kpar[2] ){ 
        	rad = 1.8957e-22;
        }
    	else if ( temp > kpar[1] ){
        	rad = 8.8669e-17/temp;
        }
    	else if ( temp > kpar[0] ){
        	rad = 1.0909e-31*pow(temp,2);
        }
    	else if ( temp >= kpar[0] ){ 
        	rad = 1.0909e-31*pow(temp,2);
        }
    	else{
        	rad = 0.0;
        }
	}
	else
	{
		//RTV loss function
		if (temp > kpar[5]){ 
        	rad = pow(10.,-17.73)/pow(temp,0.667);
        }
    	else if (temp > kpar[4] ){
        	rad = pow(10.,-21.94);
        }
    	else if (temp > kpar[3] ){
        	rad = pow(10.,-10.4)/pow(temp,2);
        }
    	else if (temp > kpar[2] ){ 
        	rad = pow(10.,-21.2);
        }
    	else if ( temp > kpar[1] ){
        	rad = pow(10.,-31.0)*pow(temp,2);
        }
    	else if ( temp >= kpar[0] ){ 
        	rad = pow(10.,-21.85);
        }
    	else{
        	rad = 0.0;
        }
	}
	
	return rad;
}



/***********************************************************************************

FUNCTION NAME: ebtel_calc_tr_dem

FUNCTION_DESCRIPTION: This function calculates the differential emission measure in the transition region using the new method described in Paper II. The method used here is to rewrite the energy equation for the transition region in terms of dtds such thatit is quadratic in dtds. The parameters aaa,bbb,ccc are the coefficients on dtds in this quadratic equation. 

INPUTS:
	tdem--temperature in the TR for which we calculate the DEM
	n--density at time t_i
	v--velocity at time t_i
	p--pressure at time t_i
	L--loop half length
	sc--scale height at temperature T_i
	rad_dem--radiative loss in the transition region
	c_a--array holding coefficients root_c2, c3, c4
	f_a--array holding f,f_eq,cf
	option--option to either use old (1) or new method (0)
	
OUTPUTS:
	dem_tr--differential emission measure in the transition region

***********************************************************************************/

double ebtel_calc_tr_dem(double tdem, double n, double v, double p, double L, double sc, double rad_dem, double c_a[3], double f_a[3], int option)
{
	double dem_tr;	//Declare what will be returned by the function

	//First check to see what method we are using to calculate the DEM in the TR
	if(option==0)
	{
		/*********New Method*********/
		//Declare variables
		double a;
		double b;
		double c;
		double p2kt2;
		double dtds1;
		double dtds2;
		double dtds;
	
		//Calculate necessary coefficients for quadratic formula
		a = KAPPA_0*pow(tdem,1.5);
		b = -5.0*K_B*n*v;
		p2kt2 = pow((p*exp(2.0*sin(PI/5.0)*L/(PI*sc))/(2*K_B*tdem)),2);	//(p/2kT)^2; calculate this here for convenience
		c = -p2kt2*rad_dem;
		dtds1 = (-b + sqrt(pow(b,2) - 4.0*a*c))/(2.0*a);
		dtds2 = (-b - sqrt(pow(b,2) - 4.0*a*c))/(2.0*a);
		dtds = ebtel_max_val(dtds1,dtds2);
		dem_tr = 2.0*p2kt2/dtds;		//factor of 2 for both legs of the loop
	}
	else
	{
		/*********Old Method*********/
		//Declare variables
		double dem_ev;
		double dem_con;
		double dem_eq;
		double root_tdem = sqrt(tdem);
		double root_rad_dem = sqrt(rad_dem);
		double fourth_tdem = pow(tdem,0.25);
		
		//Approximation to TR reg. DEM when evaporation dominates
		dem_ev = (c_a[0]/n)*(c_a[0]*pow(p,2)/root_tdem)/v;
		//Approximation to TR reg. DEM when condensation dominates
		dem_con = c_a[1]*n*v/rad_dem;
		//Approximation to TR reg. DEM under equilibrium conditions
		dem_eq = c_a[2]*p/(root_rad_dem*fourth_tdem);

		//Calculate DEM for the transition region
		dem_tr = 2*(f_a[0]*dem_ev + f_a[2]*dem_eq - f_a[1]*dem_con)/(f_a[0] + f_a[2] - f_a[1]); //factor of 2 for both legs
	}
	
	return dem_tr;

}

/***********************************************************************************

FUNCTION NAME: ebtel_heating

FUNCTION_DESCRIPTION: This function sets up the heating array to be used in the heatings
of our coronal loop. Currently, three different types of heating are available: triangular,
square, or Gaussian pulse. 

INPUTS:
	time--time array for our model
	tau--timestep
	h_nano--maximum nanoflare heating
	t_pulse_half--half the duration of the heating pulse
	heating_shape--indicates which heating profile will be used
	
OUTPUTS:
	heat--heating array

***********************************************************************************/

double * ebtel_heating(double time[], double tau, double h_nano, double t_pulse_half, int n, int heating_shape)
{
	//Declare variables
	double h_back;
	double h_thick;
	double t_start;
	double t_pulse;
	double t_end;
	double t_mid;
	double t_m;
	double t_h;
	double n_start;
	double n_end;
	double n_mid;
	double *heat = malloc(sizeof(double[n]));
	int i;

	//First set some general parameters
	h_back = 3.4e-6;
	h_thick = 0;
	t_start = 0;
	t_pulse = 2*t_pulse_half;
	t_mid = t_start + t_pulse_half;
	t_end = t_start + t_pulse;
	
	//Set indices
	n_start = t_start/tau;
	n_end = n_start + t_pulse/tau;
	n_mid = n_start + t_pulse/tau/2;
	
	//Choose which heating model to use
	//1--triangular pulse (recommended, used in Paper I,II)
	//2--square pulse 
	//3--Gaussian pulse
	
	if(heating_shape == 1)
	{
		for(i=0;i<n;i++)
		{
			//Triangular Pulse
			if(i <= n_start)
			{
				heat[i] = h_back;
			}
			else if(i > n_start && i <= n_mid)
			{
				heat[i] = h_back + h_nano*((time[i] - (t_start +tau))/t_mid);
			}
			else if(i > n_mid && i <= n_end )
			{
				heat[i] = h_back - h_nano*(time[i] - t_end)/t_mid;
			}
			else
			{
				heat[i] = h_back;
			}
		}
    }
	else if(heating_shape == 2)
	{
		for(i=0;i<n;i++)
		//Square Pulse
		if(time[i] <= (t_start + tau))
		{
			heat[i] = h_back;
		}
		else if(time[i] > t_start && time[i] < (t_end + tau))
		{
			heat[i] = h_back + h_nano;
		}
		else
		{
			heat[i] = h_back;
		}
    }
	else
	{
		//Gaussian
		//set some parameters especially for the Gaussian heating
		h_back = 3e-5;
		t_m = 2000;
		t_h = 40;
		h_nano = 1;
		for(i=0;i<n;i++)
		{		
			heat[i] = h_back + h_nano*exp(-pow((time[i] - t_m),2)/(2*pow(t_h,2)));
		}
	}

	//Return the heat pointer
	return heat;
}
 
