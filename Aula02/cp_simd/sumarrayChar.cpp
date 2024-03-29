
#include <stdio.h>
#include <time.h>
#include <cstdlib>


#define SIZE   (1000)
#define REPEAT (100000)

/**
 * sumarray using mmx instructions 
 * */
void sumarray_mmx( char *a, char *b, char *c, int size )
{

  for (int i=0;i<size;i+=8) {
  // for (int i=0;i<size;i++) {
    __asm__ volatile
        ( // instruction         comment          
        "\n\t movq     %1,%%mm0     \t#"
        "\n\t movq     %2,%%mm1     \t#"
        "\n\t paddb    %%mm0,%%mm1    \t#"
        "\n\t movq     %%mm1,%0     \t#"
        : "=m" (c[i])      // %0
        : "m"  (a[i]),     // %1 
          "m"  (b[i])      // %2
        );  
  }

   __asm__("emms" : : );
}

/**
 * sumarray using using the movdqa and paddd instructions and SSE registers
 * */
void sumarray_sse( char *a, char *b, char *c, int size )
{

  // for (int i=0;i<size;i+=4) {
  for (int i=0;i<size;i+=16) {
    __asm__ volatile
        ( // instruction         comment          
        "\n\t movq     %1,%%xmm0     \t#"
        "\n\t movq    %2,%%xmm1     \t#"
        "\n\t paddb    %%xmm0,%%xmm1    \t#"
        "\n\t movq    %%xmm1,%0     \t#"
        : "=m" (c[i])      // %0
        : "m"  (a[i]),     // %1 
          "m"  (b[i])      // %2
        ); 
  }
   __asm__("emms" : : );
}

/**
 * sumarray using classic code 
 * */
void sumarray( char *a, char *b, char *c, int size )
{
  for (int i=0;i<size;i++) {
      c[i]=a[i]+b[i];
  }
}

/**
 * print array
 * */
void print_array(char *a, int size)
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
void initArrays( char *a, char *b, char *c, int size )
{
    // for (int i=0; i< SIZE; i++) {
    //     a[i]=(i<<16)+1;
    //     b[i]=0xffff;
    //     c[i]=0;
    // }
    for (int i=0; i< SIZE; i++) {
        a[i]=(i<<4)+1;
        b[i]=0xf;
        c[i]=0;
    }
}


/**
 * test summation functions
 */
int main(int argc, char* argv[])
{
    int size = atoi(argv[1]);
    printf("size: %d\n",size);

    char a[size];
    char b[size],c[size];

    int n, nelemsum;

    clock_t init, end;

    //initialize arrays
    nelemsum=size;
    initArrays(a,b,c,nelemsum);

    // test classic code
    init = clock();
    for(n=0;n<REPEAT;n++)
        sumarray(a,b,c,nelemsum);
    end = clock();

    print_array(c,12);

    float sum_time =  (end-init)/(CLOCKS_PER_SEC*1.0);

    printf("sumarray time = %f\n",sum_time);

    //initialize arrays
    initArrays(a,b,c,nelemsum);

    // test mmx code
    init = clock();
    for(n=0;n<REPEAT;n++)
        sumarray_mmx(a,b,c,nelemsum);
    end = clock();

    print_array(c,12);

    float sumMMX_time =  (end-init)/(CLOCKS_PER_SEC*1.0);

    printf("sumarray time = %f\n", sumMMX_time);
    printf("speedup = %.2fx\n", sum_time/sumMMX_time);

    printf("\n");

    // test sse code
    init = clock();
    for(n=0;n<REPEAT;n++)
        sumarray_sse(a,b,c,nelemsum);
    end = clock();

    print_array(c,12);

    float sumSSE_time =  (end-init)/(CLOCKS_PER_SEC*1.0);

    printf("sumarray time = %f\n", sumSSE_time);
    printf("speedup = %.2fx\n", sum_time/sumSSE_time);

    return 0;
}

    
