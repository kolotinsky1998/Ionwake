CC=icc

CFLAGS=-c -fast

all: project

project: main.o  kinetic.o converter.o poisson.o output.o
	$(CC) main.o kinetic.o converter.o poisson.o output.o -o project
main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp
kinetic.o: kinetic.cpp
	$(CC) $(CFLAGS) kinetic.cpp
converter.o: converter.cpp
	$(CC) $(CFLAGS) converter.cpp
poisson.o: poisson.cpp
	$(CC) $(CFLAGS) poisson.cpp
output.o: output.cpp
	$(CC) $(CFLAGS) output.cpp


clean:
	rm -rf *.o project
