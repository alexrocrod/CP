#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <pthread.h>

void *printch(void *ach){

    char ch;
    ch = (char) ach;
    printf("%c\n",ch);
    pthread_exit(NULL);
}

int main( int argc, char** argv) 
{
   if(argc !=2) return 1;

    int sizestr = strlen(argv[1]);
    // pthread_t ths[sizestr];

    pthread_t ths[1024];
    int rc;

    for (int i = 0; i < sizestr; i++)
    {
        rc = pthread_create(&ths[i], NULL, printch, (void *) argv[1][i]);

        if(rc){
            printf("Error: unable to create thread, %d\n",i);
            exit(1);
        }
    }
    for (int i = 0; i < sizestr; i++)
    {
        pthread_join(ths[i],NULL);
    }

    return 0;   
}