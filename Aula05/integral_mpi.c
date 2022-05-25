#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main (int argc, char *argv[]) {

	int i, n;
	double PI25DT = 3.141592653589793238462643;
	double sum, h, pi, x;
	
	int numprocs, myid;
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	while (1) {
		if (myid ==0){
			printf("\nIntroduza número de intervalos: \n");
			scanf(" %d", &n);
		}
	
		MPI_Bcast(&n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
	
		double mypi=0;				
		if (n == 0) {
			break;
		} else {
			h = 1./(double)n;
			sum = 0.0;
			
			for (i = myid; i < n; i += numprocs) 
			{
				x = h*((double)i+0.5);
				sum += 4./(1+x*x);
			}

			mypi = sum*h;
			MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD);
					
			if (myid ==0){
				printf("        O valor correto de pi é %.25f\n",PI25DT);
				printf("O valor de pi é aproximadamente %.25f (n=%d intervalos).\nO erro da aproximação é %.2e%%.\n", pi, n, fabs(PI25DT-pi)/PI25DT*100);
			}
		}	
	}

	MPI_Finalize();
	
	return 0;
}
