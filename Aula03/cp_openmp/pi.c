/*
adapted from "A “Hands-on” Introduction to OpenMP" Tim Mattson Intel Corp. timothy.g.mattson@intel.com 

This program will numerically compute the integral of

                  4/(1+x*x) 
				  
from 0 to 1.  The value of this integral is pi -- which 
is great since it gives us an easy way to check the answer.

The is the original sequential program.  It uses the timer
from the OpenMP runtime library

History: Written by Tim Mattson, 11/99.

*/


#include <stdio.h>
#include <omp.h>

static long num_steps = 1000000000;
double step;

int main ()
{
	int i;
	double x, pi, sum = 0.0;
	double start_time, run_time;

	step = 1.0/(double) num_steps;

			
	start_time = omp_get_wtime();

	for (i=0; i< num_steps; i++){
		x = (i+0.5)*step;
		sum = sum + 4.0/(1.0+x*x);
	}

	pi = step * sum;
	run_time = omp_get_wtime() - start_time;
	printf("\n pi with %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time);

// -----------------------------------------------------------------------------------------------
	
	i = 0;
	x = 0.0;
	pi = 0.0;
	sum = 0.0;
	
	double run_time_OMP;

	step = 1.0/(double) num_steps;
			
	start_time = omp_get_wtime();

	#pragma omp parallel //num_threads(10)
	{
		int rank = omp_get_thread_num();
		int my_n = num_steps / omp_get_num_threads();
		double sum_th = 0.0;
		for (i= rank*my_n; i< (rank+1)*my_n; i++){
			x = (i+0.5)*step;
			sum_th = sum_th + 4.0/(1.0+x*x);
		}
		#pragma omp atomic
		sum += sum_th;
	}
	

	pi = step * sum;
	run_time_OMP = omp_get_wtime() - start_time;
	printf("\n pi with %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time_OMP);
	printf("SPEEDUP: %f \n ",run_time/run_time_OMP);

	// -----------------------------------------------------------------------------------------------
	
	i = 0;
	x = 0.0;
	pi = 0.0;
	sum = 0.0;

	step = 1.0/(double) num_steps;
			
	start_time = omp_get_wtime();

	#pragma omp parallel //num_threads(10)
	{
		double sum_th = 0.0;
		#pragma omp for
		for (i=0; i< num_steps; i++){
			x = (i+0.5)*step;
			sum_th = sum_th + 4.0/(1.0+x*x);
		}
		#pragma omp atomic
		sum += sum_th;
	}
	

	pi = step * sum;
	run_time_OMP = omp_get_wtime() - start_time;
	printf("\n pi with %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time_OMP);
	printf("SPEEDUP: %f \n ",run_time/run_time_OMP);

	// -----------------------------------------------------------------------------------------------
	
	i = 0;
	x = 0.0;
	pi = 0.0;
	sum = 0.0;

	step = 1.0/(double) num_steps;
			
	start_time = omp_get_wtime();

	#pragma omp parallel for reduction(+:sum)
	for (i=0; i< num_steps; i++){
		x = (i+0.5)*step;
		sum = sum + 4.0/(1.0+x*x);
	}

	pi = step * sum;
	run_time_OMP = omp_get_wtime() - start_time;
	printf("\n pi with %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time_OMP);
	printf("SPEEDUP: %f \n ",run_time/run_time_OMP);
}	  

