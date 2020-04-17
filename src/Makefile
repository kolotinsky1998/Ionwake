CC=g++

CFLAGS=-Wall -O3

LIBPATH=

PFLAGS=-fopenmp

all: project

project: main.o TScheme.o
	$(CC) $(CFLAGS) $(PFLAGS) main.o TScheme.o -o $(PROJECTNAME)
main.o: main.cpp
	$(CC) $(CFLAGS) $(PFLAGS) -c main.cpp $(LIBPATH)
scheme.o: TScheme.cpp
	$(CC) $(CFLAGS) $(PFLAGS) -c TScheme.cpp $(LIBPATH)
clean:
	rm *.o
	rm $(PROJECTNAME)
