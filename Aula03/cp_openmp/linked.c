// adapted from "A “Hands-on” Introduction to OpenMP" Tim Mattson Intel Corp. timothy.g.mattson@intel.com 

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#ifndef N
#define N 5
#endif
#ifndef FS
#define FS 38
#endif

struct node {
   int data;
   int fibdata;
   struct node* next;
};

int fib(int n) {
   int x, y;
   if (n < 2) {
      return (n);
   } else {
      x = fib(n - 1);
      y = fib(n - 2);
	  return (x + y);
   }
}

void processwork(struct node* p) 
{
   int n;
   n = p->data;
   p->fibdata = fib(n);
}

struct node* init_list(struct node* p) {
   int i;
   struct node* head = NULL;
   struct node* temp = NULL;
   
   head = malloc(sizeof(struct node));
   p = head;
   p->data = FS;
   p->fibdata = 0;
   for (i=0; i< N; i++) {
      temp  =  malloc(sizeof(struct node));
      p->next = temp;
      p = temp;
      p->data = FS + i + 1;
      p->fibdata = i+1;
   }
   p->next = NULL;
   return head;
}

int main(int argc, char *argv[]) {
   double start, end;
   struct node *p=NULL;
   struct node *temp=NULL;
   struct node *head=NULL;
   
   printf("Process linked list\n");
   printf("  Each linked list node will be processed by function 'processwork()'\n");
   printf("  Each ll node will compute %d fibonacci numbers beginning with %d\n\n",N,FS);      

   p = init_list(p);
   head = p;

   start = omp_get_wtime();
   {
      while (p != NULL) {
         // #pragma omp task
         processwork(p);
         p = p->next;
      }
   }


   end = omp_get_wtime();
   p = head;
   while (p != NULL) {
      printf("%d : %d, ",p->data, p->fibdata);
      temp = p->next;
      free (p);
      p = temp;
   }  
   free (p);

   double Tsingle = end-start;
   printf("\nCompute Time: %f seconds\n\n", Tsingle);

   ////////////////////////////////////////////
   p=NULL;
   temp=NULL;
   head=NULL;  

   p = init_list(p);
   head = p;

   start = omp_get_wtime();

   struct node* ps[N+1];
   int cps = 0;
   {
      while (p != NULL) {
         ps[cps] = p;
         p = p->next;
         cps++;
      }
   }

   #pragma omp parallel firstprivate(ps)
   {
      #pragma omp for
      for(int i=0; i<N+1; i++){
         processwork(ps[i]);
      }
   }
      
   end = omp_get_wtime();
   p = head;
   while (p != NULL) {
      printf("%d : %d, ",p->data, p->fibdata);
      temp = p->next;
      free (p);
      p = temp;
   }  
   free (p);

   double Twotask = end - start;
   printf("\nCompute Time w/o task: %f seconds\n", Twotask);
   printf("Speedup: %.2fx\n\n", Tsingle/Twotask);


   ////////////////////////////////////////////
   p=NULL;
   temp=NULL;
   head=NULL;  

   p = init_list(p);
   head = p;

   start = omp_get_wtime();
   #pragma omp parallel
   {
      #pragma omp single private(p)
      {
         p = head;
         while (p != NULL) {
            #pragma omp task
            processwork(p);
            p = p->next;
         }
      }
   }

   end = omp_get_wtime();
   p = head;
   while (p != NULL) {
      printf("%d : %d, ",p->data, p->fibdata);
      temp = p->next;
      free (p);
      p = temp;
   }  
   free (p);

   Twotask = end - start;
   printf("\nCompute Time with task: %f seconds\n", Twotask);
   printf("Speedup: %.2fx\n", Tsingle/Twotask);

   return 0;
}

