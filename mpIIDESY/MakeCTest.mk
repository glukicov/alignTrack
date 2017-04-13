#make file for Mptest2.cpp (translation to C++ from mptest2.f90)

PROGNAME = C_test
SOURCES =   
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLIBS     := $(shell root-config --nonew --libs)

CPP = g++
CPPFLAGS = -std=c++0x
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

$(PROGNAME) : $(OBJECTS) $(PROGNAME).o
	$(CPP) -o $@ $(OBJECTS) $(PROGNAME).o $(LDFLAGS) $(LIBS)

%.o : %.cpp
	$(CPP) $(CPPFLAGS) -o $@ -c $<

test:
	@echo $(ROOTCFLAGS)

clean :
	-rm -f ${PROGNAME} ${OBJECTS}