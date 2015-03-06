/***********************************************************************************

FILENAME: ebtel_functions_loop.c

AUTHOR Will Barnes

DATE: created: 7 March 2014

***********************************************************************************/

//Include appropriate header file
#include "ebtel_functions.h"

/***********************************************************************************

FUNCTION NAME: ebtel_loop_solver

FUNCTION DESCRIPTION: Given some heating input, this function solves the 0-D hydrodynamic
equations spatially averaged over a loop half-length, computing the temperature, density, 
and pressure as functions of time. Given appropriate initial conditions (see below), the 
differential emission measure of the transition region and corona can also be calculated. 

INPUTS:
	ntot--total number of steps to be taken in the time integration
	loop_length--half-length of the coronal loop (cm) (top of chrmosphere to apex)
   	opt--keyword structure; see ebtel_functions.h for a complete list of all structure members

OUTPUTS:
   	param_setter--pointer structure of arrays for all loop parameters; see ebtel_functions.h for a complete list of all structure members

***********************************************************************************/

struct ebtel_params_st *ebtel_loop_solver( int ntot, double loop_length, struct Option *opt)
{
	/***********************************************************************************
								Variable Declarations
	***********************************************************************************/
	
	/*****Variable declarations*****/
	//Int
	int nk;
	int i;	//index over ntot
	int j;	//index over 451
	int k;	//index used for averaging over temporal DEM dimensions
	int flag_dem_tr;
	int j_min;
	int j_max;
	int mem_lim = ntot;
	int new_mem_lim;
	int count_reallocate;
	
	//Pointers
	double *kptr;
	double *flux_ptr;
	double *state_ptr;
	double *log_tdem_ptr;
	double *ic_ptr;
	double *dem_cor_minus;
	double *dem_tr_minus;
	double *dem_tot_minus;
	
	//Double
	double r1;
	double r1_tr;
	double r2;
	double r3;
	double r4;
	double rad;
	double sc;
	double f;
	double f_eq;
	double r12;
	double r12_tr;
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
	double time = 0;	//initialize time to zero
	
	//Array 
	double f_array[3];
	double state[3];
	double log_tdem[opt->index_dem];
	double tdem[opt->index_dem];
	double root_tdem[opt->index_dem];
	double fourth_tdem[opt->index_dem];
	double rad_dem[opt->index_dem];
	double root_rad_dem[opt->index_dem];
	
	//Two-dimensional array (dynamically allocate memory to avoid segmentation fault)
	double **dem_tr;
	dem_tr = (double **)malloc(ntot*sizeof(double *));
	double **dem_cor;
	dem_cor = (double **)malloc(ntot*sizeof(double *));
	for(i = 0; i<ntot; i++)
	{
		dem_tr[i] = (double *)malloc(opt->index_dem*sizeof(dem_tr));
		dem_cor[i] = (double *)malloc(opt->index_dem*sizeof(dem_cor));
	}
	
	//struct
	struct rk_params par;
	struct ebtel_rka_st *adapt;
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
	param_setter->rad = malloc(sizeof(double[ntot]));
	param_setter->tau = malloc(sizeof(double[ntot]));
	
	if(strcmp(opt->usage_option,"rad_ratio") == 0 || strcmp(opt->usage_option,"dem") == 0)
	{
		if(strcmp(opt->usage_option,"rad_ratio") == 0)
		{	
			param_setter->f_ratio = malloc(sizeof(double[ntot]));
			param_setter->rad_ratio = malloc(sizeof(double[ntot]));
		}
		param_setter->logtdem = malloc(sizeof(double[opt->index_dem]));
		param_setter->dem_tr_log10mean = malloc(sizeof(double[opt->index_dem]));
		param_setter->dem_cor_log10mean = malloc(sizeof(double[opt->index_dem]));
		param_setter->dem_tot_log10mean = malloc(sizeof(double[opt->index_dem]));
	}
	
	/***********************************************************************************
									Initial Parameters
	***********************************************************************************/	
	
	//Calculate initial r values. In Klimchuk et al. (2008), these are denoted by c1,c2,c3. 
	//See ebtel_main.c documentation for correspondence between variables.
	r1 = ebtel_calc_c3();	//ratio of base to apex temperature; c3 in Klimchuk et al. (2008), Paper I from here on.
	r1_tr = r1;				//ratio of temperature at top of TR to apex temperature, can be set to 0.5 as well.
	r2 = ebtel_calc_c2();	//ratio of average to apex temperature; c2 in Paper I
	r3 = 2.0;				//ratio of TR to coronal radiative losses; c1 in Paper I
	r4 = 1.0;				//ratio of average to base velocity
	
	//Set temperature bins.
	//Lengths of the kpar array are different depending on the loss function we use.
	if (strcmp(opt->rad_option,"rk")== 0)
	{nk = 7;
	}
	else
	{nk = 6;
	}
	//Declare kpar array which will set the bins based on our loss function choice
	double kpar[nk];
	//Call the ebtel_kpar_set function and pass the pointer kptr
	kptr = ebtel_kpar_set(opt->rad_option);
	//Set the kpar array;
	for(i = 0; i < nk; i++)
	{
		kpar[i] = *(kptr + i);
	}
	
	/***********************************************************************************
						Set up DEM in Transition Region
	***********************************************************************************/
	if (strcmp(opt->usage_option,"dem") == 0 || strcmp(opt->usage_option,"rad_ratio") == 0)
	{	
		//Define temperature arrays for plotting and calculating DEM
		log_tdem_ptr = ebtel_linspace(0,opt->index_dem-1,opt->index_dem);		//set the linspace pointer using the ebtel_linspace function
		for(j = 0; j<opt->index_dem; j++)
		{
			log_tdem[j] = *(log_tdem_ptr + j)/100 + 4;	//log of T in the region in which DEM is defined
			param_setter->logtdem[j] = log_tdem[j];		//Save logtdem to our parameter structure
			tdem[j] = pow(10,log_tdem[j]);				//T in the region in which DEM is defined
			root_tdem[j] = sqrt(tdem[j]);				//T^1/2 in the region in which DEM is defined
			fourth_tdem[j] =  pow(tdem[j],0.25);			//T^1/4 in the region in which DEM is defined
		}
		
		//These coefficients will be used in the old method of calculating DEM (global variables)
		ROOT_C2 = pow((KAPPA_0/(20*K_B)),0.5)/K_B;		//Calculate the root of c2 to avoid overflow in our calculation of dem_ev
		C3 = -5*K_B;
		C4 = pow((KAPPA_0/14),0.5)/K_B;
		
		//Radiation in the transition region. This loop just calculates the radiative loss function in the TR
		for(j = 0; j<opt->index_dem; j++)
		{
			rad = ebtel_rad_loss(tdem[j],kpar,opt->rad_option);		//Calculate the radiative loss function for temperature tdem[i]
			rad_dem[j] = rad;								//Set radiative loss function in the TR
			if (tdem[j] < 1e+4)
			{
				rad_dem[j] = 1;								//Check to see if we are outside the allowed temperature range
			}
			root_rad_dem[j] = sqrt(rad_dem[j]);
		}
	}
	
	/***********************************************************************************
							Initial Static Equilibrium
	***********************************************************************************/
	/*
	2 possible methods: (A) use EBTEL equilibrium (recommended) (B) use scaling laws (set by mode in opt structure)
	*/
	
	//Calculate initial temperature, density, pressure, and velocity.
	ic_ptr = ebtel_calc_ic(kpar,r3,loop_length,opt);
	r3 = *(ic_ptr + 0);
	rad = *(ic_ptr + 1);
	t = *(ic_ptr + 2);
	n = *(ic_ptr + 3);
	p = *(ic_ptr + 4);
	v = *(ic_ptr + 5);
	
	//Free the pointer used in the static equilibrium calculation
	free(ic_ptr);
	ic_ptr = NULL;
	
	//Set remaining initial parameters before iterating in time
	ta = t/r2;
	sc = ebtel_calc_lambda(t);
	na = n*r2*exp(-2.0*loop_length/(PI*sc)*(1.0-sin(PI/5.0)));
	pa = 2*K_B*na*ta;
	
	//Display initial conditions when using initial conditions from static equilibrium or scaling laws
	if(strcmp(opt->ic_mode,"st_eq") == 0 || strcmp(opt->ic_mode,"scaling") == 0)
	{
		printf("************************************************************************************\n");
		printf("            		Initial Conditions		                       \n");
		printf("************************************************************************************\n");
		printf("T(t = 0) = %f K\n",t);
		printf("n(t = 0) = %f cm^-3\n",n);
		printf("p(t = 0) = %f dyne cm^-2\n",p);
		printf("r3(t = 0) = %f\n",r3);
		printf("\n");
	}
		
	//Set the initial values of our parameter structure
	param_setter->temp[0] = t;
	param_setter->ndens[0] = n;
	param_setter->press[0] = p;
	param_setter->vel[0] = v;
	param_setter->tapex[0] = ta;
	param_setter->napex[0] = na;
	param_setter->papex[0] = pa;
	param_setter->rad[0] = rad;
	param_setter->coeff_1[0] = r3;
	param_setter->heat[0] = ebtel_heating(time,opt);
	param_setter->time[0] = time;
	param_setter->tau[0] = opt->tau;
	
	/***********************************************************************************
							Time-dependent Heating
	***********************************************************************************/
	
	//Set the structure members of the par structure. This is used in both the RK and Euler methods.
	par.L = loop_length;
	par.kpar = kptr;
	par.r12 = r1/r2;
	par.r2 = r2;
	
	//Set the initial timestep from opt structure.
	//This will be static if we are not using our adaptive solver
	tau = opt->tau;
	
	//Initialize the counter
	i = 0;
	
	//Begin the loop over the timesteps
	while(time < opt->total_time)
	{
		//Update the parameter structure
		par.q1 = ebtel_heating(time,opt);
		par.q2 = ebtel_heating(time+tau,opt);
		
		//Update time
		time = time + tau;
		
		//Save time and heat to main data structure
		param_setter->heat[i+1] = par.q2;
		param_setter->time[i+1] = time;
		
		//Set up non-thermal electron flux if the nt_ebeam option was selected
		//Extensive testing for this option has not been carried out with EBTEL-C
		if(strcmp(opt->usage_option,"nt_ebeam")==0)
		{
			par.flux_nt = loop_length*par.q2;
		}
		else
		{
			par.flux_nt = 0.;
		}

		//Calculate radiative loss
		rad = ebtel_rad_loss(t,kpar,opt->rad_option);
		param_setter->rad[i+1] = rad;

		//Calculate coefficients r1, r2, r3 (c3, c2, c1)
		r3 = ebtel_calc_c1(t,n,loop_length,rad);
		par.r3 = r3;
		param_setter->coeff_1[i+1] = r3;
		r2 = ebtel_calc_c2();
		r1 = ebtel_calc_c3();
		r12 = r1/r2;
		r12_tr = r1_tr/r2;

		//Calculate thermal conduction
		flux_ptr = ebtel_calc_thermal_conduction(t, n, loop_length, rad, r3, opt->sat_limit, opt->heat_flux_option);
		f = *(flux_ptr + 0);
		f_eq = *(flux_ptr + 1);
		free(flux_ptr);
		flux_ptr = NULL;
		
		//Set thermal conduction and equilibrium heat flux in the par structure
		par.f = f;
		par.f_eq = f_eq;
		
		//Calculate pv quantity to be used in velocity calculation (enthalpy flux)
		pv = 0.4*(f_eq - f - par.flux_nt);
		
		/*****Step parameters forward in time using (1) RK method or (0) Euler stepper method*****/
		
		//Update the state vector
		state[0] = p;
		state[1] = n;
		state[2] = t;
		
		if(strcmp(opt->solver,"euler")==0)	//Euler solver
		{	
			//Call the Euler routine
			state_ptr = ebtel_euler(state,tau,par,opt);
		}
		else if(strcmp(opt->solver,"rk4")==0)	//RK routine
		{	
			//Call the RK routine
			state_ptr = ebtel_rk(state,3,param_setter->time[i],tau,par,opt);	
		}
		else if(strcmp(opt->solver,"rka4")==0)
		{	
			//Call the adaptive RK routine
			adapt = ebtel_rk_adapt(state,3,param_setter->time[i],tau,par,opt);
			//Set the state vectore and timestep
			state_ptr = adapt->state;
			tau = adapt->tau;
		}
		else
		{
			printf("Invalid solver option.\n");
			exit(0);
		}

		//Update p,n,t,tau and save to structure
		p = *(state_ptr + 0);
		param_setter->press[i+1] = p;
		n = *(state_ptr + 1);
		param_setter->ndens[i+1] = n;
		t = *(state_ptr + 2);
		param_setter->temp[i+1] = t;
		param_setter->tau[i+1] = tau;
		
		//Free memory used by the state pointer. Free the adapt structure if we are using the adapt method.
		if(strcmp(opt->solver,"rka4")==0)
		{
			//Free the structure as it will be malloc'd on the next go around
			free(adapt->state);
			adapt->state = NULL;
			free(adapt);
		}
		else
		{
			//Free the state ptr as it will be malloc'd again on the next go around
			free(state_ptr);
			state_ptr = NULL;
		}
		
		//calculate new velocity
		v = pv/p; 			
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
		
		/***********************************************************************************
							Differential Emission Measure (DEM) Calculation
		***********************************************************************************/
		
		//Check usage variable to determine whether we are calculating TR DEM
		if(strcmp(opt->usage_option,"dem") == 0 || strcmp(opt->usage_option,"rad_ratio") == 0)
		{	
			//Transition region
			if(r12_tr*t > tdem[opt->index_dem-1])
			{
				printf(" Transition region T = %e K outside of DEM range\n",r12_tr*t);
				exit(0);
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
			
			for(j=0; j<opt->index_dem; j++)
			{
				//Check to see whether we are in the TR. If so, calculate dem_TR. Note: r12_tr*t[i] = T_0
				if( tdem[j] < r12_tr*t )
				{
					//Make call to function that calculates DEM for TR
					dem_tr[i][j] = ebtel_calc_tr_dem(tdem[j],n,v,p,loop_length,sc,rad_dem[j],f_array,opt->dem_option);
					
					//Check whether the dem is less than zero and set the flag if the condition holds
					if(dem_tr[i][j] < 0)
					{
						flag_dem_tr = 1;
						printf("***** Negative DEM at i = %d\n", i);
						break;
					}
					
					//Check whether to calculate the flux ratio
					if(strcmp(opt->usage_option,"rad_ratio")==0)
					{
						//Check whether it is classical or dynamic to set the f_ratio value
						if(strcmp(opt->heat_flux_option,"classical") == 0)
						{
							//Ratio of classical flux to equilibrium flux
							f_ratio = f/f_eq;
						}
						else if(strcmp(opt->heat_flux_option,"limited")==0)
						{	
							//Ratio of classical flux to saturated flux
							f_ratio = (-TWO_SEVENTHS*KAPPA_0*pow(t/r2,SEVEN_HALVES)/loop_length)/(1.0*-1.5*pow(K_B,1.5)/pow(M_EL,0.5)*n*pow(t,1.5));
						}
						else
						{
							printf("Invalid heat flux calculation option.\n");
							exit(0);
						}
					}
				}
				else
				{
					//If outisde of the TR, set the TR DEM to zero
					dem_tr[i][j] = 0.0;
				}
			}
			
			//If we tripped the flag_dem_tr, then we need to reset dem_tr
			if(flag_dem_tr == 1)
			{
				for(j=0; j<opt->index_dem; j++)
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
			t_max = ebtel_max_val(t/r2,1.0e+4);
			t_min = ebtel_max_val(t*(2.0 - 1/r2),1.0e+4);
			
			j_max = (log10(t_max) - 4.0)*100;
			j_min = (log10(t_min) - 4.0)*100;
			
			//Calculate total emission measure
			em = 2*pow(n,2)*loop_length;			//factor of 2 for both legs
			
			delta_t = 1.0e+4*(pow(10.0,((j_max + 0.5)/100.0)) - pow(10.0,((j_min - 0.5)/100.0)));
        	dem0 = em/delta_t;
        	
        	//Set every value between j_min and j_max DEM in the corona to dem0
        	for(j=0; j<opt->index_dem; j++)
        	{
				if(j <= j_max && j >= j_min)
				{
	        		dem_cor[i][j] = dem0;
				}
				else
				{
					dem_cor[i][j] = 0.0;
				}
        	}
        	
			//Check whether to calculation the radiation ratio
        	//Transition region radiation losses based on DEM
        	if(strcmp(opt->usage_option,"rad_ratio") == 0)
        	{
        		rad_loss = 0;		//Initialize rad_loss to 0
        		
        		//Sum up radiative losses in the transition region
        		for(j=0; j<opt->index_dem; j++)
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
		i++;
		
		//Check if we need to reallocate memory
		if(i == mem_lim-1 && time < opt->total_time-tau)
		{
			//Tell the user that memory is being reallocated
			printf("Reached current memory limit.Reallocating...\n");
			
			//Increment the reallocation counter
			count_reallocate = count_reallocate + 1;
			//Update the memory limit
			new_mem_lim = 2*mem_lim;
			
			//Call the reallocation function for the param_setter structure
			ebtel_reallocate_mem(mem_lim,new_mem_lim,param_setter,opt);
			//Call the reallocation function for the two-dimensional arrays
			dem_tr = ebtel_reallocate_two_d_array(dem_tr,mem_lim,new_mem_lim,opt->index_dem);
			dem_cor = ebtel_reallocate_two_d_array(dem_cor,mem_lim,new_mem_lim,opt->index_dem);
			
			//Tell the user the new memory size and the number of reallocations performed
			printf("The new memory limit is %d\n",new_mem_lim);
			printf("Number of memory reallocations: %d\n",count_reallocate);
			//Update the memory limit
			mem_lim = new_mem_lim;
		}
	}
	
	//End of loop
	
	/***********************************************************************************
						Format Data and Free Pointers
	***********************************************************************************/
	
	//Set structure field to be used in file writer
	param_setter->i_max = i;
	
	//Format DEM data to be printed to file.
	//Take weighted time average for each temperature bin
	if(strcmp(opt->usage_option,"dem") == 0 || strcmp(opt->usage_option,"rad_ratio") == 0)
	{
		for(j = 0; j < opt->index_dem; j++)
		{
			//Set the last entry that was left empty by the loop on t
			dem_tr[param_setter->i_max-1][j] = dem_tr[param_setter->i_max-2][j];
			dem_cor[param_setter->i_max-1][j] = dem_cor[param_setter->i_max-2][j];
			
			//Malloc reduced dimension pointers
			dem_cor_minus = malloc(sizeof(double)*param_setter->i_max);
			dem_tr_minus = malloc(sizeof(double)*param_setter->i_max);
			dem_tot_minus = malloc(sizeof(double)*param_setter->i_max);
			
			//Create a single time array for each slice in j
			for(k = 0; k<param_setter->i_max; k++)
			{
				dem_cor_minus[k] = dem_cor[k][j];
				dem_tr_minus[k] = dem_tr[k][j];
				dem_tot_minus[k] = dem_tr[k][j] + dem_cor[k][j];
			}
			
			//Compute the time average at each temperature bin and take the log10
			param_setter->dem_cor_log10mean[j] = log10(ebtel_weighted_avg_val(dem_cor_minus,param_setter->i_max,param_setter->tau));
			param_setter->dem_tr_log10mean[j] = log10(ebtel_weighted_avg_val(dem_tr_minus,param_setter->i_max,param_setter->tau));
			param_setter->dem_tot_log10mean[j] = log10(ebtel_weighted_avg_val(dem_tot_minus,param_setter->i_max,param_setter->tau));
			
			//Free the reduced dimension pointers; they get malloc'd on the next iteration
			free(dem_cor_minus);
			dem_cor_minus = NULL;
			free(dem_tr_minus);
			dem_tr_minus = NULL;
			free(dem_tot_minus);
			dem_tot_minus = NULL;
			
			//Make sure that we have no negative numbers as a result of log10(0.0); -infinity *should* be ignored when plotting
			if(param_setter->dem_cor_log10mean[j] < 0.0)
			{
				param_setter->dem_cor_log10mean[j] = -INFINITY;
			}
			if(param_setter->dem_tr_log10mean[j] < 0.0)
			{
				param_setter->dem_tr_log10mean[j] = -INFINITY;
			}
			if(param_setter->dem_tot_log10mean[j] < 0.0)
			{
				param_setter->dem_tot_log10mean[j] = -INFINITY;
			}
		}
	}
	
	//Free up memory used by ebtel_kpar_set and ebtel_linspace functions
	free(kptr);
	kptr = NULL;
	if(strcmp(opt->usage_option,"dem")==0 || strcmp(opt->usage_option,"rad_ratio")==0)
	{
		free(log_tdem_ptr);
		log_tdem_ptr = NULL;
	}
	for(i=0;i<mem_lim;i++)
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
	
	//Exit and return the structure that has been set appropriately
	return param_setter;
}

/***********************************************************************************

FUNCTION NAME: ebtel_kpar_set

FUNCTION DESCRIPTION: This function sets the kpar array to be used later in calculations
of the radiative loss function. These are essentially the temperature bins used in 
calculating the radiative loss function. Either the RK or RTV method is used.

INPUTS:
	rad_option--use either the Raymond-Klimchuk or Rosner-Tucker-Vaina method 
OUTPUTS:
	*kpar--pointer for the kpar array

***********************************************************************************/

double * ebtel_kpar_set(char *rad_option)
{	
	//Check option input to decide which method to use
	if (strcmp(rad_option,"rk") == 0)
	{	
		double *kpar = malloc(sizeof(double[7]));
		//Raymond-Klimchuk Loss function
		kpar[0] = 1.0e+4;
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
		kpar[0] = pow(10.0,4.3);
        kpar[1] = pow(10.0,4.6);
        kpar[2] = pow(10.0,4.9);
        kpar[3] = pow(10.0,5.4);
        kpar[4] = pow(10.0,5.75);
        kpar[5] = pow(10.0,6.3);
        
        return kpar;
	}
}

/***********************************************************************************

FUNCTION NAME: ebtel_rad_loss

FUNCTION_DESCRIPTION: This function calculates the radiative loss function using either the Raymond-Klimchuk (RK) method or the Rosner-Tucker-Vaiana (RTV) method. Adopted from pro radLoss from the original EBTEL IDL code.

INPUTS:
	temp--temperature (K)
	kpar--holds the kpar structure if it has been previously set; placeholder otherwise.
	rad_option--option to use the Rosner-Tucker-Vaiana (rtv) or the Raymond-Klimchuk (rk) method.
	
OUTPUTS:
	rad--radiative loss function

***********************************************************************************/

double ebtel_rad_loss( double temp, double kpar[], char *rad_option)
{
	//Declare rad 
	double rad;
	
	//Find out which method is being used
	if (strcmp(rad_option,"rk") == 0)
	{
		//RK loss function
    	if ( temp > kpar[6] ){ 
        	rad = 1.96e-27*pow(temp,0.5);
        }
    	else if ( temp > kpar[5] ){ 
        	rad = 5.4883e-16/temp;
        }
    	else if ( temp > kpar[4] ){
        	rad = 3.4629e-25*(pow(temp,1./3.));
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
    	else if ( temp >= kpar[0] ){
        	rad = 1.0909e-31*pow(temp,2.);
        }
    	else{
        	rad = 0.0;
        }
	}
	else if(strcmp(rad_option,"rtv") == 0)
	{
		//RTV loss function
		if (temp > kpar[5]){ 
        	rad = pow(10.,-17.73)/pow(temp,2./3.);
        }
    	else if (temp > kpar[4] ){
        	rad = pow(10.,-21.94);
        }
    	else if (temp > kpar[3] ){
        	rad = pow(10.,-10.4)/pow(temp,2.);
        }
    	else if (temp > kpar[2] ){ 
        	rad = pow(10.,-21.2);
        }
    	else if ( temp > kpar[1] ){
        	rad = pow(10.,-31.0)*pow(temp,2.);
        }
    	else if ( temp >= kpar[0] ){ 
        	rad = pow(10.,-21.85);
        }
    	else{
        	rad = 0.0;
        }
	}
	else
	{
		printf("Invalid radiative loss function option.\n");
		exit(0);
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
	f_a--array holding f,f_eq,cf
	option--option to either use old or new method for calculating the TR DEM
	
OUTPUTS:
	dem_tr--differential emission measure in the transition region

***********************************************************************************/

double ebtel_calc_tr_dem(double tdem, double n, double v, double p, double L, double sc, double rad_dem, double f_a[3], char *option)
{
	double dem_tr;	//Declare what will be returned by the function

	//First check to see what method we are using to calculate the DEM in the TR
	if(strcmp(option,"new")==0)
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
		p2kt2 = pow((p*exp(2.0*sin(PI/5.0)*L/(PI*sc))/(2.0*K_B*tdem)),2.0);	//(p/2kT)^2; calculate this here for convenience
		c = -p2kt2*rad_dem;
		dtds1 = (-b + sqrt(pow(b,2.0) - 4.0*a*c))/(2.0*a);
		dtds2 = (-b - sqrt(pow(b,2.0) - 4.0*a*c))/(2.0*a);
		dtds = ebtel_max_val(dtds1,dtds2);
		dem_tr = 2.0*p2kt2/dtds;		//factor of 2 for both legs of the loop
	}
	else if(strcmp(option,"old")==0)
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
		dem_ev = (ROOT_C2/n)*(ROOT_C2*pow(p,2)/root_tdem)/v;
		//Approximation to TR reg. DEM when condensation dominates
		dem_con = C3*n*v/rad_dem;
		//Approximation to TR reg. DEM under equilibrium conditions
		dem_eq = C4*p/(root_rad_dem*fourth_tdem);

		//Calculate DEM for the transition region
		dem_tr = 2.*(f_a[0]*dem_ev + f_a[2]*dem_eq - f_a[1]*dem_con)/(f_a[0] + f_a[2] - f_a[1]); //factor of 2 for both legs
	}
	else
	{
		printf("Invalid TR DEM calculation option\n");
		exit(0);
	}
	
	return dem_tr;

}
