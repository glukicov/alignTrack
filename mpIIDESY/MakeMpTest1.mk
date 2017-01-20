#make file for Mptest2.cpp (translation to C++ from mptest2.f90)

PROGNAME      = MpTest1
SOURCES       = mptest1_port.cpp 
INCLUDES      = 
OBJECTS       = $(patsubst %.cpp, %.o, $(SOURCES))
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLIBS     := $(shell root-config --nonew --libs)
CFLAGS       += $(ROOTCFLAGS)
LIBS         += $(ROOTLIBS)
#
LIBS         += -lMinuit
LDFLAGS       = -O

$(PROGNAME):    $(OBJECTS)
		g++ -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS) -std=c++11

%.o : %.cpp $(INCLUDES)
	g++ ${CFLAGS} -c  -g -o $@ $< -std=c++11 -I . 

test:
	@echo $(ROOTCFLAGS)

clean:	
	-rm -f ${PROGNAME} ${OBJECTS}
