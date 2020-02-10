CC=g++

CFLAGS=-Wall -O3

LIBPATH=

PFLAGS=-fopenmp

all: project

project: main.o  scheme.o converter.o
	$(CC) $(CFLAGS) $(PFLAGS) main.o scheme.o converter.o -o $(PROJECTNAME) 
main.o: main.cpp
	$(CC) $(CFLAGS) $(PFLAGS) -c main.cpp $(LIBPATH)
scheme.o: scheme.cpp
	$(CC) $(CFLAGS) $(PFLAGS) -c scheme.cpp $(LIBPATH)
converter.o: converter.cpp
	$(CC) $(CFLAGS) -c converter.cpp
clean:
	rm *.o
	rm $(PROJECTNAME)
