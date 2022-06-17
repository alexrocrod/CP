
.PHONY: all clean

TARGETS = testDiffs harrisDetectorCudaTexture

all: $(TARGETS)

testDiffs: testDiffs.c
	g++ -O3 -ICommon $^ -o $@

harrisDetectorCudaTexture: harrisDetectorCudaTexture.cu
	nvcc -arch=sm_30 -O3 -ICommon $^ -o $@

clean:
	rm -f $(TARGETS) 
	rm -f referenceCudaTexture.pgm resultCudaTexture.pgm