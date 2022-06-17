
#include <stdio.h>
#include <time.h>
#include <cstdlib>

#define SIZE   (1000*4)
#define REPEAT (100000)

/**
 * sumarray using mmx instructions 
 * */
// int sumarray_mmx( int *a, int size )
// {
//   int res;
//   for (int i=0;i<size;i+=2) {
//     __asm__ volatile
//         ( // instruction         comment          
//         "\n\t movq     %1,%%mm0     \t#"
//         "\n\t movq     %2,%%mm1     \t#"
//         // "\n\t paddd    %%mm0,%%mm1    \t#"
//         "\n\t paddusb    %%mm0,%%mm1    \t#"
//         "\n\t movq     %%mm1,%0     \t#"
//         : "=m" (res)      // %0
//         );  
//   }

//   __asm__("emms" : : );

//   return res;
// }

/**
 * sumarray using using the movdqa and paddd instructions and SSE registers
 * */
int sumarray_sse( int *a, int size )
{
  int res = 0, temp = 0;
  for (int i=0;i<size;i+=4) {
    __asm__ volatile
        ( // instruction         comment          
        "\n\t movdqa     %1,%%xmm0     \t#"
        "\n\t haddps    %%xmm0,%%xmm1    \t#"
        "\n\t movdqa     %%xmm1,%0     \t#"
        : "=m" (temp)    // 0
        : "m" (a[i])    // 1
        ); 
    
    res += temp;
  }

  __asm__("emms" : : );

  return res;
}

/**
 * sumarray using classic code 
 * */
int sumarray( int *a, int size )
{
  int res = 0;
  for (int i=0;i<size;i++) {
      res += a[i];
  }
  return res; 
}

/**
 * print array
 * */
void print_array(int *a, int size)
{
    printf("base10: ");
    for (int i=0; i < size; i++) {
    printf("%10d",a[i]);
    }
    printf("\nbase16: ");
    for (int i=0; i < size; i++) {
    printf("%10x",a[i]);
    }
    printf("\n");
}

/**
 * init arrays
 * */
void initArrays( int *a, int size, int supSize )
{
    for (int i=0; i< size; i++) {
        a[i]= i; //(i<<16)+1;
    }
    for (int i=size; i< supSize; i++) {
        a[i]= 0; //(i<<16)+1;
    }
}


/**
 * test summation functions
 */
int main(int argc, char* argv[])
{
    int size = atoi(argv[1])*4;
    printf("size: %d\n",size);
    // int size = SIZE;
    int a[size];
    
    int n, nelemsum, res;

    clock_t init, end;

    //initialize arrays
    nelemsum=size;
    int supSize = nelemsum + (4-nelemsum%4);

    printf("s = %d, SS = %d\n",nelemsum,supSize);

    initArrays(a,nelemsum, supSize);

    // test classic code
    init = clock();
    for(n=0;n<REPEAT;n++)
        res = sumarray(a,nelemsum);
    end = clock();

    printf("res = %d\n",res);

    float sum_time =  (end-init)/(CLOCKS_PER_SEC*1.0);

    printf("sumarray time = %f\n\n",sum_time);

    // //initialize arrays
    // initArrays(a,nelemsum);

    // // test mmx code
    // init = clock();
    // for(n=0;n<REPEAT;n++)
    //     res = sumarray_mmx(a,nelemsum);
    // end = clock();

    // printf("res = %d\n",res);

    // float sumMMX_time =  (end-init)/(CLOCKS_PER_SEC*1.0);

    // printf("sumarray time = %f\n", sumMMX_time);
    // printf("speedup = %.2fx\n\n", sum_time/sumMMX_time);

    //initialize arrays
    initArrays(a, nelemsum, supSize);

    res = 0;

    // test sse code
    init = clock();
    for(n=0;n<REPEAT;n++)
        res = sumarray_sse(a,supSize);
    end = clock();

    printf("res = %d\n",res);

    float sumSSE_time =  (end-init)/(CLOCKS_PER_SEC*1.0);

    printf("sumarray time = %f\n", sumSSE_time);
    printf("speedup = %.2fx\n", sum_time/sumSSE_time);

    return 0;
}

    
