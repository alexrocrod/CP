#include <stdio.h>
#include <math.h>

int main (int argc, char *argv[]) {

	int i, n;
	double PI25DT = 3.141592653589793238462643;
	double sum, h, pi, x;

	while (1) {
		printf("\nIntroduza número de intervalos: ");
		scanf(" %d", &n);
			
		if (n == 0) {
			break;
		} else {
			h = 1./(double)n;
			sum = 0.;

			for (i = 0; i < n; i += 1) {
				x = h*((double)i+0.5);
				sum += 4./(1+x*x);
			}
			pi = sum*h;
			printf("O valor de pi é aproximadamente %f (n=%d intervalos).\nO erro da aproximação é %.2e%%.\n", pi, n, fabs(PI25DT-pi)/PI25DT*100);
		}	
	}
	return 0;
}
