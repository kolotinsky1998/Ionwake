CC=mpicxx

CFLAGS=-Wall -O3

LIBPATH=

PFLAGS=-fopenmp -lm

all: project

#PROJECTNAME=test

project: main.o TScheme.o
	$(CC) $(CFLAGS) $(PFLAGS) main.o TScheme.o -o $(PROJECTNAME)
main.o: main.cpp
	$(CC) $(CFLAGS) $(PFLAGS) -c main.cpp $(LIBPATH)
TScheme.o: TScheme.cpp
	$(CC) $(CFLAGS) $(PFLAGS) -c TScheme.cpp $(LIBPATH)
clean:
	rm *.o
	rm $(PROJECTNAME)
