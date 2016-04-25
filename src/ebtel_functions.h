/************************************************************************************

FILENAME: ebtel_functions.h

AUTHOR Will Barnes

DATE: created: 7 March 2014

DESCRIPTION: This file gives the function prototypes for functions in defined in ebtel_functions.c. Global variables are also declared here. See ebtel_functions.c and ebtel_main.c for more information on EBTEL and its original implementation in the IDL programming language (see Klimchuk et al., 2008). 

************************************************************************************/

#ifndef EBTEL_FUNCTIONS_H
#define EBTEL_FUNCTIONS_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

//Declare global variables
double K_B;
double KAPPA_0;
double M_P;
double M_EL;
double PI;
double TWO_SEVENTHS;
double SEVEN_HALVES;
double TWO_THIRDS;
double ROOT_C2;
double C3;
double C4;
int NK;
double KPAR[6];

//Declare structures
struct Option {
	int total_time;
	int index_dem;
	int num_events;
	int sample_rate;
	double loop_length;
	double energy_nt;
	double T0;
	double n0;
	double sat_limit;
	double h_nano;
	double h_back;
	double t_pulse_half;
	double t_start;
	double mean_t_start;
	double std_t_start;
	double alpha;
	double amp0;
	double amp1;
	double tau;
	double rka_error;
	double r3_rad_0, r3_eqm_0;
	double *amp;
	double *t_start_array;
	double *t_end_array;
	char *heat_flux_option;
	char *dem_option; 
	char *rad_option; 
	char *usage_option; 
	char *solver;
	char *ic_mode;
	char *print_plasma_params;
	char *output_file;
	char *heating_shape;
	char *t_start_switch;
	char *t_end_switch;
	char *amp_switch;
	char *r3_loss_correction;
	char *r3_grav_correction;
};
struct ebtel_params_st {
		int i_max;
		double *heat;
		double *time;
		double *tau;
		double *temp;
		double *ndens;
		double *press;
		double *vel;
		double *tapex;
		double *napex;
		double *papex;
		double *coeff_1;
		double *logtdem;
		double *cond;
		double *rad_cor;
		double *rad;
		double *f_ratio;
		double *rad_ratio;
		double *dem_tr_log10mean;
		double *dem_cor_log10mean;
		double *em_cor_log10mean; 
		double *dem_tot_log10mean;
};
struct rk_params {
		double L;
		double r12;
		double r2;
		double r3;
		double q1;
		double q2;
		double flux_nt;
		double f;
		double f_eq;
};
struct ebtel_rka_st {
		double tau;
		double *state;
};
struct box_muller_st{
	double z;
	double z_save;
	int flag;
};

//Declare prototype for ebtel_loop_solver of type struct *
struct ebtel_params_st *ebtel_loop_solver( int, double, struct Option *);

//Declare prototype for ebtel_kpar_set of type void
double * ebtel_kpar_set(char *);

//Declare prototype for ebtel_rad_loss of type double
double ebtel_rad_loss( double, char *);

//Declare prototype for ebtel_calc_thermal_conduction of type double *
double * ebtel_calc_thermal_conduction(double, double, double, double, double, double, char *);

//Declare prototype for ebtel_calc_c1 of type double
double ebtel_calc_c1(double, double, double, double, struct Option *);

//Declare prototype for ebtel_calc_c2 of type double
double ebtel_calc_c2(void);

//Declare prototype for ebtel_calc_c3 of type double
double ebtel_calc_c3(void);

//Declare prototype for ebtel_calc_lambda of type double
double ebtel_calc_lambda(double);

//Declare prototype for ebtel_linspace of type void
double * ebtel_linspace( int, int, int);

//Declare prototype for ebtel_calc_tr_dem of type double
double ebtel_calc_tr_dem( double, double, double, double, double, double, double, double[], char *);

//Declare prototype for ebtel_avg_val of type double
double ebtel_avg_val(double[], int );

//Declare protoype for ebtel_max_val of type double
double ebtel_max_val(double, double);

//Declare prototype for ebtel_mem_free of type void
void ebtel_free_mem(struct ebtel_params_st *, struct Option *);

//Declare prototype for ebtel_rk of type double
double * ebtel_rk(double[], int, double, double, struct rk_params, struct Option *);

//Declare prototype for ebtel_heating of type double
double ebtel_heating(double, struct Option *);

//Declare prototype for ebtel_print_header of type void
void ebtel_print_header(int, struct Option *);

//Declare prototype for ebtel_derivs function of type double
double * ebtel_derivs(double[], double, struct rk_params, struct Option *);

//Declare prototype for ebtel_data_writer of type void
void ebtel_file_writer(struct Option *, struct ebtel_params_st *);

//Declare prototype for ebtel_rk_adapt of type struct ebtel_rka_st
struct ebtel_rka_st *ebtel_rk_adapt(double[], int, double, double, struct rk_params, struct Option *);

//Declare prototype for ebtel_min_val of type double
double ebtel_min_val(double, double);

//Declare prototype for ebtel_calc_abundance of type void
void ebtel_calc_abundance(void);

//Declare prototype for ebtel_static_eq of type double
double * ebtel_calc_ic(double, double, struct Option *);

//Declare prototype for ebtel_colon_operator of type double *
double * ebtel_colon_operator(double, double, double);

//Declare prototype for ebtel_weighted_avg_val of type double
double ebtel_weighted_avg_val(double *, int, double *);

//Declare prototype for ebtel_box_muller of type struct box_muller_st
struct box_muller_st *ebtel_box_muller(double,double,double,int);

//Declare prototype for ebtel_rand_limit of type double
double ebtel_rand_limit(double);

//Declare prototype for ebtel_power_law of type double
double ebtel_power_law(double,double,double,double);

//Declare prototype for ebtel_bubble_sort of type double *
double * ebtel_bubble_sort(double[],int);

//Declare prototype for ebtel_heating_profiles of type double
double ebtel_heating_profile(double,double,double,double,struct Option *);

//Declare prototype for ebtel_heating_config of type void
void ebtel_heating_config(struct Option *, char *);

//Declare prototype for ebtel_count_events of type int
int ebtel_count_events(struct ebtel_params_st *, struct Option *);

//Declare prototype for ebtel_input_setter of type struct Option
struct Option *ebtel_input_setter(char *);

//Declare prototype for ebtel_xml_reader of type char *
char *ebtel_xml_reader(xmlNodePtr,char *, char *);

//Declare prototype for ebtel_reallocate_mem of type void
void ebtel_reallocate_mem(int, int, struct ebtel_params_st *, struct Option *);

//Declare prototype for ebtel_reallocate_two_d_array of type double **
double **ebtel_reallocate_two_d_array(double **, int, int, int);

#endif