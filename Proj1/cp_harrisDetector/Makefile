
.PHONY: all clean

TARGETS = harrisDetectorCuda harrisDetectorOpenMP testDiffs harrisDetectorCudaConst harrisDetectorCudaTexture

all: $(TARGETS)

harrisDetectorCuda: harrisDetectorCuda.cu
	nvcc -arch=sm_30 -O3 -ICommon $^ -o $@

harrisDetectorOpenMP: harrisDetectorOpenMP.c
	g++ -fopenmp -O3 -ICommon $^ -o $@

testDiffs: testDiffs.c
	g++ -O3 -ICommon $^ -o $@

harrisDetectorCudaConst: harrisDetectorCudaConst.cu
	nvcc -arch=sm_30 -O3 -ICommon $^ -o $@

harrisDetectorCudaTexture: harrisDetectorCudaTexture.cu
	nvcc -arch=sm_30 -O3 -ICommon $^ -o $@

clean:
	rm -f $(TARGETS) 
	rm -f referenceCuda.pgm resultCuda.pgm referenceOpenMP.pgm resultOpenMP.pgm referenceCudaTexture.pgm resultCudaTexture.pgm referenceCudaConst.pgm resultCudaConst.pgm
