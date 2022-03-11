#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <pthread.h>

void *dot1(void *axs){

    float *xs;
    xs = (float*) axs;
    // printf("x1 = %f, x2 = %f\n",xs[0], xs[1]);
    xs[2] = xs[0]*xs[1];
    pthread_exit(NULL);
}

int main( int argc, char** argv) 
{
    float vec1[] = {2.5, 10};
    float vec2[] = {10, 1};
    printf("x1 = %f, y1 = %f\n",vec1[0], vec1[1]);
    printf("x2 = %f, y2 = %f\n",vec2[0], vec2[1]);

    float res = 0.0f;
    pthread_t ths[2];

    float xs[3];
    xs[0] = vec1[0];
    xs[1] = vec2[0];
    xs[2] = 0.0f;

    float ys[3];
    ys[0] = vec1[1];
    ys[1] = vec2[1];
    ys[2] = 0.0f;

    int rc = pthread_create(&ths[0], NULL, dot1, (void *)xs);
    if(rc){
        printf("Error: unable to create thread 0\n");
        exit(1);
    }

    rc = pthread_create(&ths[1], NULL, dot1, (void *)ys);
    if(rc){
        printf("Error: unable to create thread 1\n");
        exit(1);
    }

    for (int i = 0; i < 2; i++){
        pthread_join(ths[i],NULL);
    }

    res = xs[2] + ys[2];
    printf("res = %f\n",res);
    
    return 0;   
}