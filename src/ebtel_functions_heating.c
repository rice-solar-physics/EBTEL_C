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
 
FUNCTION NAME: ebtel_heating_config

FUNCTION DESCRIPTION: This function configures the start times, end times, and amplitudes
of the heating events as specified by the user.

INPUTS:
	opt--pointer to option structure

OUTPUTS:

***********************************************************************************/

void ebtel_heating_config(struct Option *opt, char *filename)
{
	//Variable declarations
	int bm_flag = 0;
	int i;
	
	double x1,x2;
	double tmp,save;
	double limit = 1.;
	double *sort_ptr1;
	double *sort_ptr2;
	
	char temp[64];
	
	//Declare doc variable
	xmlDocPtr doc;
	//Create the document tree
	doc = xmlParseFile(filename);
	//Create a node for the root of the tree
	xmlNodePtr root = doc->children;
	
	//Declare structure for Box-Muller algorithm
	struct box_muller_st *bm_st;
	
	//Declare amplitude and start time arrays
	double amp[opt->num_events];
	double t_start_array[opt->num_events];
	double t_end_array[opt->num_events];

	//Reserve memory for amplitude and start time arrays in opt structure
	opt->t_start_array = malloc(sizeof(double[opt->num_events]));
	opt->amp = malloc(sizeof(double[opt->num_events]));
	opt->t_end_array = malloc(sizeof(double[opt->num_events]));

	//Seed the random number generator
	srand(time(NULL));
	
	//Calculate the start times and amplitudes
	//Begin loop to set start times
	for(i=0;i<opt->num_events;i++)
	{
		//Set random numbers for either start time or amplitudes
		if(strcmp(opt->t_start_switch,"normal") ==0 || strcmp(opt->amp_switch,"power_law") == 0)
		{
			//Initialize the two random variables
			x1 = ebtel_rand_limit(limit);
			x2 = ebtel_rand_limit(limit);
		}
	
		//Use uniformly or normally distributed start times
		if(strcmp(opt->t_start_switch,"uniform") == 0)
		{
			//Start times separated by two pulse durations (following Reep et al. 2013)
			t_start_array[i] = opt->t_start + 2.*i*(2*opt->t_pulse_half);
		}
		else if(strcmp(opt->t_start_switch,"normal") == 0)
		{
			//Use the Box-Muller method to do the normal distribution
			bm_st = ebtel_box_muller(x1,x2,save,bm_flag);
			tmp = bm_st->z;
			save = bm_st->z_save;
			bm_flag = bm_st->flag;

			//Save the 'denormalized' normally distributed start time
			t_start_array[i] = opt->std_t_start*tmp + opt->mean_t_start;
		
			//Free the structure
			free(bm_st);
			bm_st = NULL;
		}
		else if(strcmp(opt->t_start_switch,"file") == 0)
		{	
			//Set temp variable to search tree
			sprintf(temp,"start_time_array%d",i);
			//Set value of the start time array
			t_start_array[i] = atof(ebtel_xml_reader(root,temp,NULL));	
		}
		else
		{
			printf("Invalid heating start time option. Choose either uniform, normal, or file.\n");
			exit(0);
		}
	
		//Use uniform amplitudes, amplitudes given by power law distribution, or read them in from config file
		if(strcmp(opt->amp_switch,"uniform") == 0)
		{
			amp[i] = opt->h_nano;
		}
		else if(strcmp(opt->amp_switch,"power_law") == 0)
		{
			//Compute the amplitude according to a power-law distribution
			amp[i] = ebtel_power_law(opt->amp0,opt->amp1,x1,opt->alpha);
		}
		else if(strcmp(opt->amp_switch,"file") == 0)
		{
			//Set temp variable to search tree
			sprintf(temp,"amp_array%d",i);
			//Set value of the start time array
			amp[i] = atof(ebtel_xml_reader(root,temp,NULL));
		}
		else
		{
			printf("Invalid heating amplitude option. Choose either uniform, power_law, or file.\n");
			exit(0);
		}
		
		//Use uniform pulse times (as defined by configuration file) or read in pulse times from config file
		if(strcmp(opt->t_end_switch,"uniform") == 0)
		{
			//Set array of pulse times from parameter file
			t_end_array[i] = 2*opt->t_pulse_half + t_start_array[i];
		}
		else if(strcmp(opt->t_end_switch,"file") == 0)
		{
			//Set temp variable to search tree
			sprintf(temp,"end_time_array%d",i);
			//Set value of the start time array
			t_end_array[i] = atof(ebtel_xml_reader(root,temp,NULL));
		}
		else
		{
			printf("Invalid heating end time option. Choose either uniform or file.\n");
			exit(0);
		}
	}
	
	//Free the XML document tree
	xmlFreeDoc(doc);

	//If the start times are out of order, sort start and end times in ascending order.
	if(strcmp(opt->t_start_switch,"normal") == 0)
	{
		sort_ptr1 = ebtel_bubble_sort(t_start_array,opt->num_events);
		sort_ptr2 = ebtel_bubble_sort(t_end_array,opt->num_events);	
		for(i=0;i<opt->num_events;i++)
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
	for(i=0; i<opt->num_events; i++)
	{
		opt->t_start_array[i] = t_start_array[i];
		opt->t_end_array[i] = t_end_array[i];
		opt->amp[i] = amp[i];
	}
}

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
	double heat;
	double heat_temp;
		
	//Set the heating as the background heating
	heat = opt->h_back;
	
	//Check all heating intervals to see if we have fallen into one of them
	//Test all timing intervals
	for(i=0;i<opt->num_events;i++)
	{
		//Check if we are inside the heating pulse interval
		if(t >= *(opt->t_start_array + i) && t <= (*(opt->t_end_array + i) ) )
		{
			//If so, call the heating profile function to generate the correct pulse
			heat_temp = ebtel_heating_profile(t,*(opt->t_start_array + i),*(opt->t_end_array + i),*(opt->amp + i),opt);
			heat = heat + heat_temp;
		}
	}
	
	//Return the heating value
	return heat;
}

/***********************************************************************************

FUNCTION NAME: ebtel_heating_profile

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

double ebtel_heating_profile(double t, double t_start, double t_end, double h_nano, struct Option *opt)
{
	//Variable declarations and definitions
	double t_pulse = t_end - t_start;
	double t_mid = t_start + t_pulse/2.;
	double heat;
	
	//Choose which heating model to use
	//1--triangular pulse
	//2--square pulse 
	//3--Gaussian pulse
	//If you insist on hardcoding a heating function, it should be added here.
	
	if(strcmp(opt->heating_shape,"triangle") == 0)
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
	else if(strcmp(opt->heating_shape,"square") == 0)
	{
		//Square pulse
		heat = h_nano;
		
    }
	else if(strcmp(opt->heating_shape,"gaussian") == 0)
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

FUNCTION NAME: ebtel_count_events

FUNCTION_DESCRIPTION: This function counts the number of events that actually occur
within the simulation. This may differ from the number specified in the config file
due to random start times not falling inside the simulated t domain.

INPUTS:
	struct ebtel_params_st params_final--structure holding all resulting arrays
	struct Option opt--structure holding start and end time arrays
	
OUTPUTS:
	int num_events--number of heating events simulated
***********************************************************************************/

int ebtel_count_events(struct ebtel_params_st *params_final,struct Option *opt)
{
	//Declare variables
	int i;
	int N_events = 0;
	double start_time;
	double end_time;
	double start_event;
	double end_event;
	
	//Get start time and end time
	start_time = *(params_final->time + 0);
	end_time = *(params_final->time + params_final->i_max-1);
	
	//Check each event to see if it occured in the given domain
	for(i=0;i<opt->num_events;i++)
	{
		//Set event start and end times
		start_event = *(opt->t_start_array + i);
		end_event = *(opt->t_end_array + i);
		//Check the event
		if(start_event >= start_time && start_event < end_time && end_event <= end_time && end_event > start_time)
		{
			//Increment the event counter
			N_events = N_events + 1;
		}
	}
	
	//Return actual number of events
	return N_events;
}
