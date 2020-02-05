CC=icpc

CFLAGS=-ggdb -Wall -xCORE-AVX2 -ipo -O3 -no-prec-div -qopt-prefetch -ffinite-math-only -qopt-mem-layout-trans=3 

LIBPATH=-I../c_arrays-master/mdarrayXX.hpp/ -I../c_arrays-master/mdarrayXX.hpp/src-generated

PFLAGS=-qopenmp

all: project

project: main.o  scheme.o converter.o nik-log-txt-tmr-exit.o
	$(CC) $(CFLAGS) $(PFLAGS) main.o scheme.o converter.o nik-log-txt-tmr-exit.o -o $(PROJECTNAME) 
main.o: main.cpp
	$(CC) $(CFLAGS) $(PFLAGS) -c main.cpp $(LIBPATH)
scheme.o: scheme.cpp
	$(CC) $(CFLAGS) $(PFLAGS) -c scheme.cpp $(LIBPATH)
converter.o: converter.cpp
	$(CC) $(CFLAGS) -c converter.cpp
nik-log-txt-tmr-exit.o:
	$(CC) $(CFLAGS) -c nik-log-txt-tmr-exit.cpp $(LIBPATH)
clean:
	rm *.o
	rm $(PROJECTNAME)
