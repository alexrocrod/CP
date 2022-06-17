
.PHONY: all clean

TARGETS = testDiffs harrisDetectorCudaConst

all: $(TARGETS)

testDiffs: testDiffs.c
	g++ -O3 -ICommon $^ -o $@

harrisDetectorCudaConst: harrisDetectorCudaConst.cu
	nvcc -arch=sm_30 -O3 -ICommon $^ -o $@

clean:
	rm -f $(TARGETS) 
	rm -f referenceCudaConst.pgm resultCudaConst.pgm
