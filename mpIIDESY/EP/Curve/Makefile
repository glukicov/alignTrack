Target  = curve.exe
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
all:$(Target)

curve.exe: curve.C
	g++ -std=c++0x -o $@ curve.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm curve.exe

