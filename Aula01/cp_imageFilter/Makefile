
TARGETS = imageFilter, tchar

CFLAGS = -Wall -O0

all: $(TARGETS)

imageFilter: imageFilter.o
	gcc -o imageFilter imageFilter.o

tchar: tchar.o
	gcc -o tchar tchar.o -lpthread

clean:
	rm -f $(TARGETS) *.o
