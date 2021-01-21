# Directories
OBJ=/Users/zhouyan17/Desktop/bumphunter/obj
SRC=/Users/zhouyan17/Desktop/bumphunter/src
INC=/Users/zhouyan17/Desktop/bumphunter/src/inc

# Setting ROOT variables
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

LIBS += -lTMVA

# Setting compiler flags
CXX = g++ -Wall -Wformat=0 -c -g -shared -fPIC
#CXX = g++ -std=c++11 -Wall -Wformat=0 -c -g -shared -fPIC
#CXX = g++ -Wall -Wformat=0
CXXFLAGS = -I$(ROOTCFLAGS) -I/$(INC)

LD = g++
LDFLAGS	= $(ROOTLIBS) $(LIBS) -L$(INC) -L$(ROOTCFLAGS)

# List of sources and objects
sources=$(shell find $(SRC) -type f -name \*.cxx) runBumpHunter.cxx
objects=$(sources:$(SRC)/%.cxx=$(OBJ)/%.o)

# Default target
default: runBumpHunter

# Generic rule for objects
$(OBJ)/%.o: $(SRC)/%.cxx $(INC)/%.h 
	${CXX} ${CXXFLAGS} -c $< -o $@

# Rule for main object
runBumpHunter.o: runBumpHunter.cxx
	${CXX} ${CXXFLAGS} -c $< -o $@

# Rule for main program
runBumpHunter: $(objects)
	@echo "objects: $(objects)"
	$(LD) $(LDFLAGS) $^ -o $@

clean:
	rm -rf $(OBJ)/*.o runBumpHunter


