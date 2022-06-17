
TARGETS = imageFilter

CFLAGS = -Wall -O0

all: $(TARGETS)

imageFilter: imageFilter.o
	gcc -o imageFilter imageFilter.o -pthread

clean:
	rm -f $(TARGETS) *.o
