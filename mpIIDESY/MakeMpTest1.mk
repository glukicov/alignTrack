# Makefile for MpTest1 (translation to C++ from mptest1.f90)

PROGNAME = MpTest1
SOURCES =  Logger.cpp mptest1_port_detector.cpp random_buffer.cpp mptest1_port.cpp
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLIBS     := $(shell root-config --nonew --libs)

CPP = g++
CPPFLAGS = -std=c++11
CPPFLAGS += $(ROOTCFLAGS)

LDD = g++
LDFLAGS = 

LIBS += $(ROOTLIBS)
LIBS         += -lMinuit
LDFLAGS       = -O

//CFLAGS += $(shell root-config --cflags)
//LDFLAGS += -g $(shell $(ROOTSYS)/bin/root-config --ldflags)
//LIBS += $(shell root-config --glibs)

all : $(PROGNAME)

$(PROGNAME) : $(OBJECTS)
	$(CPP) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

%.o : %.cpp
	$(CPP) $(CPPFLAGS) -o $@ -c $<

test:
	@echo $(ROOTCFLAGS)

clean :
	-rm -f ${PROGNAME} ${OBJECTS}