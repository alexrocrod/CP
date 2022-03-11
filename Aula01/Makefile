
TARGETS = imageFilter, tchar, thDotP

CFLAGS = -Wall -O0

all: $(TARGETS)

imageFilter: imageFilter.o
	gcc -o imageFilter imageFilter.o

tchar: tchar.o
	gcc -o tchar tchar.o -pthread

thDotP: thDotP.o
	gcc -o thDotP thDotP.o -pthread

clean:
	rm -f $(TARGETS) *.o
