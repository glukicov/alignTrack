#make file for AlignTracker.cpp 

PROGNAME = AlignTracker
SOURCES =  Logger.cpp AlignTracker_methods.cpp random_buffer.cpp
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLIBS     := $(shell root-config --nonew --libs)

CPP = g++
#CPP=/usr/local/bin/g++-7
CPPFLAGS = -std=c++11
CPPFLAGS += $(ROOTCFLAGS)
#CPPFLAGS +=-static-libtsan
#CPPFLAGS += -Wstatic-float-init
CPPFLAGS += -fstack-protector-all
CPPFLAGS += -pedantic-errors
CPPFLAGS += -O3

LDD = g++
LDFLAGS = 

LIBS += $(ROOTLIBS)
LIBS         += -lMinuit
# LDFLAGS       = -O

CFLAGS += $(shell root-config --cflags)
LDFLAGS += -g $(shell $(ROOTSYS)/bin/root-config --ldflags)
LIBS += $(shell root-config --glibs)

all : $(PROGNAME)

$(PROGNAME) : $(OBJECTS) $(PROGNAME).o
	$(CPP) -o $@ $(OBJECTS) $(PROGNAME).o $(LDFLAGS) $(LIBS)

%.o : %.cpp %.h AlignTracker_methods.h
	$(CPP) $(CPPFLAGS) -c $<

test:
	@echo $(ROOTCFLAGS)

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CPP) $(SOURCES) -MM $^ -MF ./.depend;

include .depend

clean :
	-rm -f ${PROGNAME} ${OBJECTS} Tracker_con.txt Tracker_str.txt Tracker_data.bin