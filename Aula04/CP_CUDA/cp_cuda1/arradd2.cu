#include  <stdio.h> 
#include  <time.h> 

#define  SIZE (16) 
#define THPERBLOCK (4)

//  Kernel definition, see also section 2.1 of NVIDIA CUDA Programming Guide 
__global__  void arrAdd(float *A, float *B, float *C) 
{ 
    // TODO: determine id
    int id;

    id = blockIdx.x * blockDim.x + threadIdx.x;

    C[id] = A[id] + B[id]; 
} 

int  main(void) 
{ 
    float A[SIZE], B[SIZE], C[SIZE]; 
    float *devPtrA; 
    float *devPtrB; 
    float *devPtrC; 
    int memsize = SIZE * sizeof(float); 

    // Initialize arrays
    srand (time(NULL));
    for(int i=0; i < SIZE; i++) 
    {
        A[i]=rand() % 100;
        B[i]=rand() % 100;
    }

    cudaSetDevice(0);    // Select GPU device (can be 0 to 1)

    // Allocate device memory for A, B and C arrays
    cudaMalloc((void**)&devPtrA, memsize); 
    cudaMalloc((void**)&devPtrB, memsize); 
    cudaMalloc((void**)&devPtrC, memsize); 

    // Copy data (data to process) from host to device (from CPU to GPU)
    cudaMemcpy(devPtrA, A, memsize,  cudaMemcpyHostToDevice); 
    cudaMemcpy(devPtrB, B, memsize,  cudaMemcpyHostToDevice); 

    // Execute the Kernel
    // TO DO: Determine the number of blocks
    arrAdd <<<SIZE/THPERBLOCK, THPERBLOCK>>> (devPtrA,  devPtrB, devPtrC); // launch ??? blocks with 2 threads each

    // Copy data from device (results) back to host 
    cudaMemcpy(C, devPtrC, memsize,  cudaMemcpyDeviceToHost); 

    // Show results
    printf("     A      B       C\n");
    for (int i=0; i < SIZE; i++) 
    {
        printf("%2d: %4.1f + %4.1f = %5.1f\n", i, A[i], B[i], C[i]); 
    }

    // Free device memory
    cudaFree(devPtrA); 
    cudaFree(devPtrB); 
    cudaFree(devPtrC); 
} 

