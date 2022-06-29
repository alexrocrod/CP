
#include <stdio.h>
#include <time.h>
#include <cstdlib>

#define SIZE   (1000*4)
#define REPEAT (100000)

/**
 * sumarray using mmx instructions 
 * */
void sumarray_mmx( int *a, int *b, int *c, int size )
{

  for (int i=0;i<size;i+=2) {
    __asm__ volatile
        ( // instruction         comment          
        "\n\t movq     %1,%%mm0     \t#"
        "\n\t movq     %2,%%mm1     \t#"
        // "\n\t paddd    %%mm0,%%mm1    \t#"
        "\n\t paddusb    %%mm0,%%mm1    \t#"
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
void sumarray_sse( int *a, int *b, int *c, int size )
{

  for (int i=0;i<size;i+=4) {
    __asm__ volatile
        ( // instruction         comment          
        "\n\t movdqa     %1,%%xmm0     \t#"
        "\n\t movdqa     %2,%%xmm1     \t#"
        "\n\t paddd    %%xmm0,%%xmm1    \t#"
        "\n\t movdqa     %%xmm1,%0     \t#"
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
void sumarray( int *a, int *b, int *c, int size )
{
  for (int i=0;i<size;i++) {
      c[i]=a[i]+b[i];
  }
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
void initArrays( int *a, int *b, int *c, int size )
{
    for (int i=0; i< size; i++) {
    //     a[i]=(i<<16)+1;
    //     b[i]=0xffff;
        c[i]=0;
    }
    a[0] = 0xbbb0;
    a[1] = 0xff00;
    a[2] = 0xc678;
    a[3] = 0xa574;

    b[0] = 0x0afb;
    b[1] = 0x5106;
    b[2] = 0xa098;
    b[3] = 0x1b1f;
}


/**
 * test summation functions
 */
int main(int argc, char* argv[])
{
    // int size = atoi(argv[1])*4;
    int size = 1*4;
    // printf("size: %d\n",size);
    // // int size = SIZE;
    int a[size];
    int b[size],c[size];
    
    int n, nelemsum;

    clock_t init, end;

    //initialize arrays
    nelemsum=size;
    initArrays(a,b,c,nelemsum);

    // int a[] = {0xbbb0,0xff00,0xc678,0xa574};
    // int b[] = {0x0afb,0x5106,0xa098,0x1b1f};
    // int c[] = {0,0,0,0};

    // test classic code
    init = clock();
    for(n=0;n<REPEAT;n++)
        sumarray(a,b,c,nelemsum);
    end = clock();

    int nprint = nelemsum;
    if (nelemsum>12) nprint = 12;
    
    print_array(c,nprint);

    float sum_time =  (end-init)/(CLOCKS_PER_SEC*1.0);

    printf("sumarray time = %f\n",sum_time);

    
    //initialize arrays
    initArrays(a,b,c,nelemsum);

    // test mmx code
    init = clock();
    for(n=0;n<REPEAT;n++)
        sumarray_mmx(a,b,c,nelemsum);
    end = clock();

    
    print_array(c,nprint);

    float sumMMX_time =  (end-init)/(CLOCKS_PER_SEC*1.0);

    printf("sumarray time = %f\n", sumMMX_time);
    printf("speedup = %.2fx\n", sum_time/sumMMX_time);

    printf("\n");

    // test sse code
    init = clock();
    for(n=0;n<REPEAT;n++)
        sumarray_sse(a,b,c,nelemsum);
    end = clock();

    print_array(c,nprint);

    float sumSSE_time =  (end-init)/(CLOCKS_PER_SEC*1.0);

    printf("sumarray time = %f\n", sumSSE_time);
    printf("speedup = %.2fx\n", sum_time/sumSSE_time);

    return 0;
}

    
