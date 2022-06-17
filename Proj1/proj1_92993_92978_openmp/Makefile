
.PHONY: all clean

TARGETS = harrisDetectorOpenMP testDiffs

all: $(TARGETS)

harrisDetectorOpenMP: harrisDetectorOpenMP.c
	g++ -fopenmp -O3 -ICommon $^ -o $@

testDiffs: testDiffs.c
	g++ -O3 -ICommon $^ -o $@

clean:
	rm -f $(TARGETS) 
	rm -f referenceOpenMP.pgm resultOpenMP.pgm
