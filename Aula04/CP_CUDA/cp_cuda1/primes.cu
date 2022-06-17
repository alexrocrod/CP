#include  <stdio.h> 
#include  <time.h> 

#define  NROWS (1024) 
#define  NCOLS (8) 
#define  SIZE (NROWS*NCOLS) 


int compare(int *a1, int *a2);


//  Kernel definition, see also section 2.1 of NVIDIA CUDA Programming Guide 
__global__  void primes(int *A, int *B) 
{ 
    // TODO: determine id
    int id;

    id = (blockIdx.y * blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;


    if(id < SIZE)
    {
        if (A[id] == 0 || A[id] == 1)
        {
            B[id] = 0;
        }
        for (int i = 2; i <= A[id] / 2 && !B[id]; ++i) {
            if (A[id] % i == 0) 
            {
                B[id] = 1;
            }
        }
    }
} 


int primesHost(int A) 
{ 
    int B = 0;
    if (A == 0 || A == 1)
    {
        B = 0;
    }
    for (int i = 2; i <= A / 2 && !B; ++i) {
        if (A % i == 0) 
        {
            B = 1;
        }
    }
    return B;
} 

int  main(void) 
{ 
    int A[SIZE], D[SIZE], H[SIZE];
    int *devPtrA; 
    int *devPtrD; 
    int memsize = SIZE * sizeof(int); 
    float devExecTime;

    cudaSetDevice(0);   // Select GPU device (can be 0 to 1)

    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Initialize arrays
    srand (time(NULL));
    for(int i=0; i < SIZE; i++) 
    {
        A[i]=rand() % 100;
    }

    printf("Starting HOST...\n");

    for(int i=0; i < NCOLS; i++)
    {
        for(int j=0; j < NROWS; j++)
        {
            int id = j + i * NROWS;
            H[id] = primesHost(A[id]); /// host result...
        }
    }

    // Allocate device memory for A, B and D arrays
    cudaMalloc((void**)&devPtrA, memsize); 
    cudaMalloc((void**)&devPtrD, memsize); 

    printf("Starting DEVICE...\n");
    cudaEventRecord(start);

    // Copy data (data to process) from host to device (from CPU to GPU)
    cudaMemcpy(devPtrA, A, memsize,  cudaMemcpyHostToDevice);

    // __global__ functions are called:  Func <<< dim grid, dim block >>> (parameter); 
    dim3 dimBlock(2, 2);
    dim3 dimGrid(SIZE / dimBlock.x, SIZE / dimBlock.y);
    

    // Execute the Kernel 
    primes <<<dimGrid, dimBlock>>> (devPtrA, devPtrD); 

    // Copy data from device (results) back to host 
    cudaMemcpy(D, devPtrD, memsize,  cudaMemcpyDeviceToHost); 

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&devExecTime, start, stop); //Exec time = elapsed time
    

    // Show results
    printf("     A      B       D      H\n");
    for (int i=0; i < SIZE; i++) 
    {
        printf("%2d: %4d -> %5d [%5d]\n", i, A[i], D[i], H[i]); 
    }

    printf("\nOutput arrays (H/D) are %s\n", compare(D, H) == 1 ? "EQUAL" : "DIFFERENT");

    printf("\nDevice execution time [ms]: %7.4f\n", devExecTime);

    // Free device memory
    cudaFree(devPtrA);
    cudaFree(devPtrD); 
} 

int compare(int *a1, int *a2)
{
    int i, j, equal = 1;
    for(j=0; (j < NROWS) && equal; j++)
    {
        for(i=0; (i < NCOLS) && equal; i++)
        {
            int id = i + j * NCOLS;
            if(a1[id] != a2[id])
                equal = 0;
        }
    }
    return equal;
}


