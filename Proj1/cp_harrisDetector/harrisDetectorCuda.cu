
// Based on CUDA SDK template from NVIDIA

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>
#include <float.h>

// includes, project
#include <helper_cuda.h>
#include <helper_image.h>

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

#define MAX_BRIGHTNESS 255
 
// pixel base type
// Use int instead `unsigned char' so that we can
// store negative values.
typedef int pixel_t;


// harris detector code to run on the host
void harrisDetectorHost(const pixel_t *h_idata, const int w, const int h, 
                const int ws,               // window size
                const int threshold,        // threshold value to detect corners
                pixel_t * reference)
{
    int i,j,k,l;  // indexes in image
    int Ix, Iy;   // gradient in XX and YY
    int R;        // R metric
    int sumIx2, sumIy2, sumIxIy;

    for(i=0; i<h; i++) //height image
    {
        for(j=0; j<w; j++) //width image
        {
            reference[i*w+j]=h_idata[i*w+j]/4; // to obtain a faded background image
        }
    }

    for(i=ws+1; i<h-ws-1; i++) //height image
    {
        for(j=ws+1; j<w-ws-1; j++) //width image
        {
           sumIx2=0;sumIy2=0;sumIxIy=0;
           for(k=-ws; k<=ws; k++) //height window
              {
                  for(l=-ws; l<=ws; l++) //width window
                  {
                        Ix = ((int)h_idata[(i+k-1)*w + j+l] - (int)h_idata[(i+k+1)*w + j+l])/32;         
                        Iy = ((int)h_idata[(i+k)*w + j+l-1] - (int)h_idata[(i+k)*w + j+l+1])/32;         
                        sumIx2 += Ix*Ix;
                        sumIy2 += Iy*Iy;
                        sumIxIy += Ix*Iy;
                  }
              }

              R = sumIx2*sumIy2-sumIxIy*sumIxIy-0.05*(sumIx2+sumIy2)*(sumIx2+sumIy2);
              if(R > threshold) {
                   reference[i*w+j]=MAX_BRIGHTNESS; 
              }
        }
    }
}   

__global__ void kernel_Harris(pixel_t *dev_idata, const int w, const int h, 
                  const int ws, const int threshold, pixel_t *dev_odata)		
{
    int Ix, Iy;   // gradient in XX and YY
    int R;        // R metric
    int sumIx2=0, sumIy2=0, sumIxIy=0;

    int j = (blockIdx.x * blockDim.x) + threadIdx.x;
    int i = (blockIdx.y * blockDim.y) + threadIdx.y;

    // fade image
    dev_odata[i*w+j]=dev_idata[i*w+j]/4;

    if((i>= ws + 1) && (i < h-ws-1) && (j >= ws + 1) && (j < w-ws-1)) // only if valid interior region
    {
        for(int k=-ws; k<=ws; k++) //height window
        {
            for(int l=-ws; l<=ws; l++) //width window
            {
                Ix = ((int)dev_idata[(i+k-1)*w + j+l] - (int)dev_idata[(i+k+1)*w + j+l])/32;         
                Iy = ((int)dev_idata[(i+k)*w + j+l-1] - (int)dev_idata[(i+k)*w + j+l+1])/32;          
                sumIx2 += Ix*Ix;
                sumIy2 += Iy*Iy;
                sumIxIy += Ix*Iy;
            }
        }
    

        R = sumIx2*sumIy2-sumIxIy*sumIxIy-0.05*(sumIx2+sumIy2)*(sumIx2+sumIy2);
        
        // is a corner
        if(R > threshold) 
                dev_odata[i*w+j]=MAX_BRIGHTNESS;

    }
}

// harris detector code to run on the GPU
void harrisDetectorDevice(const pixel_t *h_idata, const int w, const int h, 
                  const int ws, const int threshold, 
                  pixel_t * h_odata)
{
    //TODO
    pixel_t *dev_idata, *dev_odata;
    pixel_t data_size;

    // full size for an image
    data_size = w * h * sizeof(pixel_t);
    printf("data_size %d\n", data_size);

    // Max number of threads
    int numThreads_x = 32;
	int numThreads_y = 32;
    
    // Resize for exact scaling
    while((w % numThreads_x ) != 0)
        numThreads_x = numThreads_x - 1;
        
	while((h % numThreads_y ) != 0)
        numThreads_y = numThreads_y - 1;

    // Number of blocks with that number of threads
    int numBlocks_x = ceil(w/numThreads_x);
	int numBlocks_y = ceil(h/numThreads_y);

    // CUDA dimensions
    dim3 dimBlock(numThreads_x, numThreads_y); //  threadsPerBlock
	dim3 dimGrid(numBlocks_x, numBlocks_y); // numBlocks

    // print dimensions
    printf ("w = %d\n", w);
	printf ("h = %d\n", h);
	printf("dimBlock = %d x %d \n", dimBlock.x, dimBlock.y);
	printf("dimGrid = %d x %d \n", dimGrid.x, dimGrid.y); 

    // memory allocation
    cudaMalloc((void **)&dev_idata, data_size);
    cudaMalloc((void **)&dev_odata, data_size);

    // copy image to device (CPU->GPU)
    cudaMemcpy(dev_idata, h_idata, data_size, cudaMemcpyHostToDevice);

    // Run corner detetion on GPU
    kernel_Harris<<<dimGrid, dimBlock>>>(dev_idata, w, h, ws, threshold, dev_odata);

    // Copy result from device to host (GPU->CPU)
    cudaMemcpy(h_odata, dev_odata, data_size, cudaMemcpyDeviceToHost);

    // free allocated memory
    cudaFree(dev_idata);
    cudaFree(dev_odata);

}

// print command line format
void usage(char *command) 
{
    printf("Usage: %s [-h] [-d device] [-i inputfile] [-o outputfile] [-r referenceFile] [-w windowsize] [-t threshold]\n",command);
}

// main
int main( int argc, char** argv) 
{

    // default command line options
    int deviceId = 0;
    char *fileIn        = (char *)"chess.pgm",
         *fileOut       = (char *)"resultCuda.pgm",
         *referenceOut  = (char *)"referenceCuda.pgm";
    unsigned int ws = 1, threshold = 500;

    // parse command line arguments
    int opt;
    while( (opt = getopt(argc,argv,"d:i:o:r:w:t:h")) !=-1)
    {
        switch(opt)
        {

            case 'd':
                if(sscanf(optarg,"%d",&deviceId)!=1)
                {
                    usage(argv[0]);
                    exit(1);
                }
                break;

            case 'i':
                if(strlen(optarg)==0)
                {
                    usage(argv[0]);
                    exit(1);
                }

                fileIn = strdup(optarg);
                break;
            case 'o':
                if(strlen(optarg)==0)
                {
                    usage(argv[0]);
                    exit(1);
                }
                fileOut = strdup(optarg);
                break;
            case 'r':
                if(strlen(optarg)==0)
                {
                    usage(argv[0]);
                    exit(1);
                }
                referenceOut = strdup(optarg);
                break;
            case 'w':
                if(strlen(optarg)==0 || sscanf(optarg,"%d",&ws)!=1)
                {
                    usage(argv[0]);
                    exit(1);
                }
                break;
            case 't':
                if(strlen(optarg)==0 || sscanf(optarg,"%d",&threshold)!=1)
                {
                    usage(argv[0]);
                    exit(1);
                }
                break;
            case 'h':
                usage(argv[0]);
                exit(0);
                break;

        }
    }

    // select cuda device
    checkCudaErrors( cudaSetDevice( deviceId ) );
    
    // create events to measure host harris detector time and device harris detector time

    cudaEvent_t startH, stopH, startD, stopD;
    checkCudaErrors( cudaEventCreate(&startH) );
    checkCudaErrors( cudaEventCreate(&stopH)  );
    checkCudaErrors( cudaEventCreate(&startD) );
    checkCudaErrors( cudaEventCreate(&stopD)  );



    // allocate host memory
    pixel_t * h_idata=NULL;
    unsigned int h,w;

    //load pgm
    if (sdkLoadPGM<pixel_t>(fileIn, &h_idata, &w, &h) != true) {
        printf("Failed to load image file: %s\n", fileIn);
        exit(1);
    }

    // allocate mem for the result on host side
    pixel_t * h_odata   = (pixel_t *) malloc( h*w*sizeof(pixel_t));
    pixel_t * reference = (pixel_t *) malloc( h*w*sizeof(pixel_t));
 
    // detect corners at host

    checkCudaErrors( cudaEventRecord( startH, 0 ) );
    harrisDetectorHost(h_idata, w, h, ws, threshold, reference);   
    checkCudaErrors( cudaEventRecord( stopH, 0 ) ); 
    checkCudaErrors( cudaEventSynchronize( stopH ) );

    // detect corners at GPU
    checkCudaErrors( cudaEventRecord( startD, 0 ) );
    harrisDetectorDevice(h_idata, w, h, ws, threshold, h_odata);   
    checkCudaErrors( cudaEventRecord( stopD, 0 ) ); 
    checkCudaErrors( cudaEventSynchronize( stopD ) );
    
    // check if kernel execution generated and error
    getLastCudaError("Kernel execution failed");

    float timeH, timeD;
    checkCudaErrors( cudaEventElapsedTime( &timeH, startH, stopH ) );
    printf( "Host processing time: %f (ms)\n", timeH);
    checkCudaErrors( cudaEventElapsedTime( &timeD, startD, stopD ) );
    printf( "Device processing time: %f (ms)\n", timeD);

    // save output images
    if (sdkSavePGM<pixel_t>(referenceOut, reference, w, h) != true) {
        printf("Failed to save image file: %s\n", referenceOut);
        exit(1);
    }
    if (sdkSavePGM<pixel_t>(fileOut, h_odata, w, h) != true) {
        printf("Failed to save image file: %s\n", fileOut);
        exit(1);
    }

    // cleanup memory
    free( h_idata);
    free( h_odata);
    free( reference);

    checkCudaErrors( cudaDeviceReset() );
}
