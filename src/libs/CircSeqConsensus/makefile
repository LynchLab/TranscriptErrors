cscSource = CS_consensus.cxx naiveRepeatSearch.cxx

BASE_LIB_DIR = ../


LIB_NAME = csc


CXX = g++
CFLAGS = -c -g


all: csc

csc: $(cscSource)
	$(CXX) $(CFLAGS) $(cscSource) -I$(BASE_LIB_DIR)
	ar rcs lib$(LIB_NAME).a $(cscSource:%.cxx=%.o)

clean:
	rm *.o *.a
