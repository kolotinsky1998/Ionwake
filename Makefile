CC := g++
#CC=icc
CFLAGS=-Ofast
#fast
PFLAGS=-fopenmp
#PFLAGS= -qopenmp

SRCDIR := src
BUILDDIR := build
SRCEXT := .cpp
TARGET := ionwake

ifeq ($(OS),Windows_NT)
	SOURCES      := $(shell FORFILES /P $(SRCDIR) /S /M *.cpp /C "CMD /C ECHO @relpath")
	SOURCES      := $(patsubst ".\\%",$(SRCDIR)\\%,$(SOURCES))
	OBJECTS      := $(SOURCES:$(SRCDIR)\\%.cpp=$(BUILDDIR)\\%.o)
	CLEANTEXT    := IF EXIST $(BUILDDIR) RMDIR /S /Q $(BUILDDIR)
	CREATEFOLDER := MD
else
	SOURCES      := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
	OBJECTS      := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
	CLEANTEXT    := rm -r $(BUILDDIR) $(TARGET)
	CREATEFOLDER := mkdir -p
endif

all: $(OBJECTS)
	@echo linking...
	$(CC) $^ -o $(TARGET) $(PFLAGS)

$(BUILDDIR)\\%.o: $(SRCDIR)\%.cpp
	@IF NOT EXIST $(dir $@) MD $(dir $@)
#	mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c -o $@ $< $(PFLAGS)

clean:
	$(CLEANTEXT)
#	rm -rf *.o project


#all: project
#
#project: main.o  kinetic.o converter.o poisson.o output.o
#	$(CC) main.o kinetic.o converter.o poisson.o output.o $(PFLAGS) -o $(PROJECTNAME)
#main.o: src/main.cpp
#	$(CC) $(CFLAGS) main.cpp $(PFLAGS)
#kinetic.o: src/kinetic.cpp
#	$(CC) $(CFLAGS) kinetic.cpp $(PFLAGS)
#converter.o: src/converter.cpp
#	$(CC) $(CFLAGS) converter.cpp
#poisson.o: src/poisson.cpp
#	$(CC) $(CFLAGS) poisson.cpp
#output.o: src/output.cpp
#	$(CC) $(CFLAGS) output.cpp
#
#
#clean:

