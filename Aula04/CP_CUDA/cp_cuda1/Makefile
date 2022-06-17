
TARGETS = arradd1 arradd2 arradd3 arradd4 primes

.PHONY: all clean

all: $(TARGETS)


arradd1: arradd1.cu
	nvcc -arch=sm_30 arradd1.cu -o arradd1

arradd2: arradd2.cu
	nvcc -arch=sm_30 arradd2.cu -o arradd2

arradd3: arradd3.cu
	nvcc -arch=sm_30 arradd3.cu -o arradd3

arradd4: arradd4.cu
	nvcc -arch=sm_30 arradd4.cu -o arradd4

primes: primes.cu
	nvcc -arch=sm_30 primes.cu -o primes

clean:
	rm $(TARGETS)
