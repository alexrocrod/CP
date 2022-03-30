// adapted from "A “Hands-on” Introduction to OpenMP" Tim Mattson Intel Corp. timothy.g.mattson@intel.com 

#include <stdio.h>
#include <omp.h>

int main ()  
{
    omp_set_num_threads(10);
    // #pragma omp parallel num_threads(10)
    #pragma omp parallel
    {
    int ID = omp_get_thread_num();

    printf("Hello (%d) ", ID);
    printf("World (%d) \n", ID);
    }
    
}
