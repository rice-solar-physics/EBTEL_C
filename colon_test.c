#include "ebtel_functions.h"

//Testing file for colon_operator function in ebtel

//Function prototype for colon_operator
double * ebtel_colon_operator(double, double, double);

int main(int argc, char *argv[])
{
	//Declare variables
	int j;
	double start;
	double end;
	double delta;
	
	if(argc > 4)
	{
		printf("Too many arguments. Exiting the program\n");
		return 1;
	}
	else
	{
		start = atof(argv[1]);
		end = atof(argv[2]);
		delta = atof(argv[3]);
	}
	
	//Calculate n
	int n = ceil((end - start)/delta);
	if(n == (end - start)/delta)
	{
		n = n+1;
	}
	printf("start = %f\n",start);
	printf("end = %f\n",end);
	printf("delta = %f\n",delta);
	printf("n = %d\n",n);
	double array1[n];
	double array2[n];
	
	//Declare the pointer
	double *array1_ptr;
	double *array2_ptr;
	
	//call the function
	array1_ptr = ebtel_colon_operator(start,end,delta);
	for(j=0;j<n;j++)
	{
		array1[j] = *(array1_ptr + j);
	}
	
	//Print some values to the screen
	printf("array1(end) = %f\n",array1[n-1]);
	printf("array1(0) = %f\n",array1[0]);
	
	//Call the linspace function
	array2_ptr = ebtel_linspace(start,end,n);
	for(j = 0; j<n; j++)
	{
		array2[j] = *(array2_ptr + j);
	}
	
	//Print some values for linspace
	printf("array2(end) = %f\n",array2[n-1]);
	printf("array2(0) = %f\n",array2[0]);
	
	//Test output
	printf("delta0 = %f\n",array1[0] - array2[0]);
	printf("deltaF = %f\n",array1[n-1] - array2[n-1]);
	
	//Free the pointer
	free(array1_ptr);
	array1_ptr = NULL;
	free(array2_ptr);
	array2_ptr = NULL;
	
	return 0;
}