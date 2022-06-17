
CFLAGS   = -Wall
CXXFLAGS = -Wall

.PHONY: all clean

all: sumarray sumarrayChar sumelems

sumarray: sumarray.o
	g++ $(CFLAGS) -o $@ $^

sumarrayChar: sumarrayChar.o
	g++ $(CFLAGS) -o $@ $^

sumelems: sumelems.o
	g++ $(CXXFLAGS) -o $@ $^

clean:
	rm -f *.o sumarray sumarrayChar sumelems

