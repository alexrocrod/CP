#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <pthread.h>

struct params{
    float *vec1, *vec2, sum;
    int start, finish;
} par1, par2;

void *dot1(void *apars){
    struct params *pars = (struct params *)apars;
    
    pars->sum = 0;
    for (int i = pars->start; i < pars->finish; i++)
    {
        pars->sum += pars->vec1[i]*pars->vec2[i];
    }
    pthread_exit(NULL);    
}

int main( int argc, char** argv) 
{
    float ola1[3] = {2.5, 10, 2};
    float ola2[3] = {10, 1, 4};

    float res = 0.0f;
    pthread_t ths[2];

    int lens = sizeof(ola1)/sizeof(float);
    par1.vec1 = ola1;
    par1.vec2 = ola2;
    par1.sum = 0;
    par1.start = 0;
    par1.finish = floor(lens/2);

    par2 = par1;
    par2.start = floor(lens/2);
    par2.finish = lens;

    printf("mid %d\n",par2.start);


    int rc = pthread_create(&ths[0], NULL, dot1, &par1);
    if(rc){
        printf("Error: unable to create thread 0\n");
        exit(1);
    }

    rc = pthread_create(&ths[1], NULL, dot1, &par2);
    if(rc){
        printf("Error: unable to create thread 1\n");
        exit(1);
    }

    pthread_join(ths[0],NULL);
    pthread_join(ths[1],NULL);

    printf("sumf1 %f\n",par1.sum);
    
    res = par1.sum + par2.sum;
    printf("res = %f\n",res);

    
    
    return 0;   
}