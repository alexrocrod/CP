
.PHONY: all clean

TARGETS = harrisDetectorCuda testDiffs
all: $(TARGETS)

harrisDetectorCuda: harrisDetectorCuda.cu
	nvcc -arch=sm_30 -O3 -ICommon $^ -o $@

testDiffs: testDiffs.c
	g++ -O3 -ICommon $^ -o $@

clean:
	rm -f $(TARGETS) 
	rm -f referenceCuda.pgm resultCuda.pgm
