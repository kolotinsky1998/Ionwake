CC=icc

CFLAGS=-c -fast

PFLAGS=

all: project

project: main.o  kinetic.o converter.o poisson.o output.o
	$(CC) $(PFLAGS) main.o kinetic.o converter.o poisson.o output.o -o project
main.o: main.cpp
	$(CC) $(CFLAGS) $(PFLAGS) main.cpp
kinetic.o: kinetic.cpp
	$(CC) $(CFLAGS) $(PFLAGS) kinetic.cpp
converter.o: converter.cpp
	$(CC) $(CFLAGS) $(PFLAGS) converter.cpp
poisson.o: poisson.cpp
	$(CC) $(CFLAGS) $(PFLAGS) poisson.cpp
output.o: output.cpp
	$(CC) $(CFLAGS) $(PFLAGS) output.cpp


clean:
	rm -rf *.o project
