#make file for Test.cpp (translation to C++ from mptest2.f90)

PROGNAME      = Test
SOURCES       = Test.cpp Logger.cpp
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
		g++ -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS) -std=c++0x

%.o : %.cpp $(INCLUDES)
	g++ ${CFLAGS} -c  -g -o $@ $< -std=c++0x -I . 

test:
	@echo $(ROOTCFLAGS)

clean:	
	-rm -f ${PROGNAME} ${OBJECTS}
