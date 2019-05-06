Target  = tmp.exe
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
all:$(Target)

tmp.exe: tmp.C
	g++ -std=c++0x -o $@ tmp.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm *exe

