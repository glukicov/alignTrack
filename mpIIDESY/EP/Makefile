C++ = g++
C_FLAGS = -g

ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)

SOURCES = $(shell find .  -maxdepth 1 ! -name "Reader.C" -name "*.C"  )

all : $(SOURCES:.C=)

clean:
	rm *.o $(SOURCES:.C=)



Reader.o: Reader.C Reader.h
	$(C++) -c $(ROOTCFLAGS) $(C_FLAGS) Reader.C -o Reader.o


# Ntuple_dict.so: LinkDef.h
# 	rootcling -f Ntuple_dict.c -p $^
# 	$(C++) -shared -c $(ROOTCFLAGS) $(C_FLAGS)  Ntuple_dict.c -o Ntuple_dict.so

# Skimmer: Ntuple_dict.so Ntuple_Classes.h Skimmer.h Skimmer.C Reader.o 
# 	$(C++) $(ROOTCFLAGS) $(C_FLAGS) Ntuple_dict.so Reader.o Ntuple_Classes.h Skimmer.C -o Skimmer $(ROOTLIBS)


##generic rule for any plotter:
%: %.C Plotter.h Reader.o
	$(C++) $(ROOTCFLAGS) $(C_FLAGS) Reader.o $<  -o $@  $(ROOTLIBS)




