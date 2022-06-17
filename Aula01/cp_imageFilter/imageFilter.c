
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <pthread.h>



#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

#define HI(num)	(((num) & 0x0000FF00) << 8)
#define LO(num)	((num) & 0x000000FF)

typedef struct 
{
    unsigned int *h_idata;
    unsigned int w;
    unsigned int h;
    float* filter;
    unsigned int fw;
    unsigned int fh; 
    unsigned int* h_odata;
    unsigned int w0;
    unsigned int h0;
} filterImageArgs;


// loads grayscale image to parray from fname, 
// allocates memory to parray and stores width, height and maxvalue of image through pwidth, pheight and pmaxgray
int loadPGM(char* fname, unsigned int** parray, unsigned int *pwidth, unsigned int *pheight, unsigned int *pmaxgray)
{
    FILE *pgmFile;
    char version[4];
    int i, j;
    int lo, hi;
    pgmFile = fopen(fname, "rb");
    if (pgmFile == NULL) {
        perror("cannot open file to read");
        return -1;
    }
    if(fgets(version, sizeof(version), pgmFile)==NULL) {
        fprintf(stderr, "error reading file\n");
        return -1;
    }
    if (strcmp(version, "P5\n")) {
        fprintf(stderr, "File format is not P5 pgm\n");
        return -1;
    }
    fscanf(pgmFile, "#%*[^\n]%*c");  // skip comment
    if(fscanf(pgmFile, "%d", pwidth)!=1) {
        fprintf(stderr, "error reading file\n");
        return -1;
    }
    if(fscanf(pgmFile, "%d", pheight)!=1) {
        fprintf(stderr, "error reading file\n");
        return -1;
    }
    if(fscanf(pgmFile, "%d", pmaxgray)!=1) {
        fprintf(stderr, "error reading file\n");
        return -1;
    }
    fgetc(pgmFile);
 
    *parray = malloc( *pwidth * *pheight * sizeof(int) );
    if (*pmaxgray > 255) {
        for (i = 0; i < *pheight; ++i) {
            for (j = 0; j < *pwidth; ++j) {
                hi = fgetc(pgmFile);
                lo = fgetc(pgmFile);
                (*parray)[i * *pwidth + j] = (hi << 8) + lo;
            }
        }
    }
    else {
        for (i = 0; i < *pheight; ++i) {
            for (j = 0; j < *pwidth; ++j) {
                lo = fgetc(pgmFile);
                (*parray)[i * *pwidth + j] = lo;
            }
        }
    }
 
    fclose(pgmFile);
    return 0;
 
}
 

// saves grayscale image on parray to fname, 
int savePGM(char* fname, unsigned int* array, unsigned int width, unsigned int height, unsigned int maxgray)
{
    FILE *pgmFile;
    int i, j;
    int hi, lo;
 
    pgmFile = fopen(fname, "wb");
    if (pgmFile == NULL) {
        perror("cannot open file to write");
        return -1;
    }
 
    fprintf(pgmFile, "P5\n");
    fprintf(pgmFile, "%d %d\n", width, height);
    fprintf(pgmFile, "%d\n", maxgray);
 
    if (maxgray > 255) {
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                hi = HI(array[i*width+j]);
                lo = LO(array[i*width+j]);
                fputc(hi, pgmFile);
                fputc(lo, pgmFile);
            }
 
        }
    }
    else {
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                lo = LO(array[i*width+j]);
                fputc(lo, pgmFile);
            }
        }
    }
    fprintf(pgmFile,"\n");
 
    fclose(pgmFile);

    return 0;
}


// loads filter coefficients from file fname, 
// allocates memory through parray and stores width and height of filter through pwidth and pheight
int loadFilter(char* fname, float** parray, unsigned int *pwidth, unsigned int *pheight)
{
    FILE* fp;

    if( (fp=fopen(fname, "r")) == NULL)
    {    
        fprintf(stderr,"Failed to open filter file %s\n",fname);
        return -1;
    }

    if(fscanf(fp,"%u %u",pwidth,pheight)!=2) {
        fprintf(stderr,"Failed to read header of filter file %s\n",fname);
        return -1;
    }

    *parray = (float *) malloc((*pwidth)*(*pheight)*sizeof(float));

    int i;
    for(i=0;i<(*pwidth)*(*pheight);i++)
    {
        if(fscanf(fp,"%f",(*parray+i))!=1) {
           fprintf(stderr,"Failed to read data of filter file %s\n",fname);
           return -1;
        }
    }

    fclose(fp);
    
    return 0;
}

// filter image in h_idata using linear filter, stores result in h_odata
void filterImage(unsigned int *h_idata, unsigned int w, unsigned int h, 
                float* filter, unsigned int fw, unsigned int fh, 
                unsigned int* h_odata, unsigned int w0, unsigned int h0)
{
    int i,j,k,l;

    int fh_2 = fh/2;
    int fw_2 = fw/2;

    h += h0;
    w += w0;

    for(i=w0; i<h; i++) //height image
    {
        for(j=h0; j<w; j++) //width image
        {
            float sum = 0;
            for(k=-fh_2; k<=fh_2; k++) //filter height
            {
                for(l=-fw_2; l<=fw_2; l++) //filter width
                {
                    if( (i+k >= 0) && (i+k < h))
                        if( (j+l >=0) && (j+l < w)) {
                            sum += h_idata[(i+k)*w + j+l]*filter[(k+fh/2)*fw + l+fw/2];         
                        }

                }
            }
            h_odata[i*w+j] = min(max(sum,0),255);
        }
    }
}   

void *filterImageTH(void *afargs)
{
    filterImageArgs *fargs = (filterImageArgs *)afargs;

    filterImage(fargs->h_idata, fargs->w, fargs->h, fargs->filter, fargs->fw, fargs->fh, fargs->h_odata, fargs->w0, fargs->h0);
    pthread_exit(NULL);
}



// print command line format
void usage(char *command) 
{
    printf("Usage: %s [-h] [-i inputfile] [-o outputfile] [-f filterfile] [-t nthreads]\n",command);
}

// main
int main( int argc, char** argv) 
{
    unsigned int* h_idata=NULL;
    unsigned int* h_odata=NULL;
    unsigned int h, w, maxgray;

    float *filter;
    unsigned int fh, fw;

    // default command line options
    char *fileIn="lena.pgm",*fileOut="lenaOut.pgm",*fileFilter="filter.txt";

    int nthreads = 2;

    // parse command line arguments
    int opt;
    while( (opt = getopt(argc,argv,"i:o:f:h:t:")) !=-1)
    {
        switch(opt)
        {
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
            case 'f':
                if(strlen(optarg)==0)
                {
                    usage(argv[0]);
                    exit(1);
                }
                fileFilter = strdup(optarg);
                break;
            case 'h':
                usage(argv[0]);
                exit(0);
                break;
            case 't':
                if(strlen(optarg)==0)
                {
                    usage(argv[0]);
                    exit(1);
                }
                nthreads = atoi(optarg);              
                printf("%d threads\n", nthreads);
                break;
        }
    }

    if (argc > optind) {
        usage(argv[0]);
        exit(0);
    }

    //load pgm
    if (loadPGM(fileIn, &h_idata, &w, &h,&maxgray)==-1) {
        printf("Failed to load image file: %s\n", fileIn);
        exit(1);
    }

    //load filter
    if(loadFilter(fileFilter, &filter, &fw, &fh)==-1)
    {
       printf("Failed to load filter file: %s\n",fileFilter);
       exit(1);
    }

    // allocate mem for the result
    h_odata   = (unsigned int*) malloc( h*w*sizeof(unsigned int));
    printf("w=%d, h=%d\n", w, h);
 
    pthread_t ths[w];
    int rc;

    struct timeval start, end;
    gettimeofday(&start, NULL);
    
    filterImage(h_idata, w, h, filter, fw, fh, h_odata,0,0); 

    gettimeofday(&end, NULL);
    double linearT =  (end.tv_sec-start.tv_sec)*1000.0 + ((double)(end.tv_usec - start.tv_usec))/1000.0;
    printf( "Linear Processing time: %f (ms)\n", linearT);
    ///------------------------------------------------------------------

    int n0ths = nthreads;

    if (nthreads > w) nthreads = w;

    // struct timeval start, end;
    gettimeofday(&start, NULL);
    {
        // Multithread filtering (vertical strips)
        unsigned int w1 = w/nthreads;
        unsigned int wlast = w-(w1*(nthreads-1));

        for (int i = 0; i < nthreads; i++){

            filterImageArgs fargs1;
            fargs1.h_idata = h_idata;
            fargs1.w = w1;
            fargs1.h = h;
            fargs1.filter = filter;
            fargs1.fw = fw;
            fargs1.fh = fh;
            fargs1.h_odata = h_odata;
            fargs1.w0 = w1*i;
            fargs1.h0 = 0;

            if (i==nthreads-1)
                fargs1.w = wlast;
                
            rc = pthread_create(&ths[i], NULL, filterImageTH, &fargs1);
            if(rc){printf("Error: unable to create thread, %d\n",i); exit(1);}

        }
        for (int i = 0; i < nthreads; i++)
        {
            pthread_join(ths[i],NULL);
        }
        
    }
    gettimeofday(&end, NULL);
    double VerticalT = (end.tv_sec-start.tv_sec)*1000.0 + ((double)(end.tv_usec - start.tv_usec))/1000.0;
    printf( "Vertical parts: Processing time: %f (ms), speedup: %.2f x\n", VerticalT, (linearT/VerticalT));

    ///------------------------------------------------------------------

    nthreads = n0ths;

    if (nthreads > h) nthreads = h;

    // struct timeval start, end;
    gettimeofday(&start, NULL);
    {
        // Multithread filtering (horizontal strips)
        unsigned int h1 = h/nthreads;
        unsigned int hlast = h-(h1*(nthreads-1));

        for (int i = 0; i < nthreads; i++){

            filterImageArgs fargs;
            fargs.h_idata = h_idata;
            fargs.w = w;
            fargs.h = h1;
            fargs.filter = filter;
            fargs.fw = fw;
            fargs.fh = fh;
            fargs.h_odata = h_odata;
            fargs.w0 = 0;
            fargs.h0 = h1*i;

            if (i==nthreads-1)
                fargs.h = hlast;

            rc = pthread_create(&ths[i], NULL, filterImageTH, &fargs);
            if(rc){printf("Error: unable to create thread, %d\n",i); exit(1);}

        }

        for (int i = 0; i < nthreads; i++)
        {
            pthread_join(ths[i],NULL);
        }
    }
    gettimeofday(&end, NULL);
    double HorizontalT = (end.tv_sec-start.tv_sec)*1000.0 + ((double)(end.tv_usec - start.tv_usec))/1000.0;
    printf( "Horizontal parts: Processing time: %f (ms), speedup: %.2f x\n", HorizontalT, (linearT/HorizontalT));

    /////---------------------------------------------------------
    nthreads = n0ths;

    if (nthreads > w) nthreads = w;

    // struct timeval start, end;
    gettimeofday(&start, NULL);
    {
        // Multithread filtering (sucesssive pixels) . by considering that successive pixels in each line 
        // will be processed by different successive threads.

        for (int j = 0; j < h; j++) // lines
        {
            for (int nth = 0; nth < w; nth += nthreads)
            {
                int maxK = min(nth + nthreads,w);
                for (int k = nth; k < maxK ; k++)
                {  
                    int i = k - nth;

                    filterImageArgs fargs;
                    fargs.h_idata = h_idata;
                    fargs.w = 1;
                    fargs.h = 1;
                    fargs.filter = filter;
                    fargs.fw = fw;
                    fargs.fh = fh;
                    fargs.h_odata = h_odata;
                    fargs.w0 = k;
                    fargs.h0 = j;

                    // printf("j=%d, k=%d, th=%d\n",j,k,i);


                    rc = pthread_create(&ths[i], NULL, filterImageTH, &fargs);
                    if(rc){printf("Error: unable to create thread, %d\n",i); exit(1);}

                    // pthread_join(ths[i],NULL);
                }
                for (int a = 0; a < maxK-nth; a++)
                {
                    pthread_join(ths[a],NULL);
                }
            }
        }
            

    }
    gettimeofday(&end, NULL);
    double ColumnsT = (end.tv_sec-start.tv_sec)*1000.0 + ((double)(end.tv_usec - start.tv_usec))/1000.0;
    printf( "Collumns: Processing time: %f (ms), speedup: %.2f x\n", ColumnsT, (linearT/ColumnsT));

    //---------------------------------------------------------
    
    // save output image
    if (savePGM(fileOut, h_odata, w, h, maxgray)==-1) {
        printf("Failed to save image file: %s\n", fileOut);
        exit(1);
    }

    // cleanup memory
    free(h_idata);
    free(h_odata);
    free(filter);
}
