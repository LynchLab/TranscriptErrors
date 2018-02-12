mapGenomeSource = circMapGenome.cxx

BASE_LIB_DIR = ../

LIB_NAME = circMapGenome


CXX = g++
CFLAGS = -c -g

all: map

map: $(mapGenomeSource)
	$(CXX) $(CFLAGS) $(mapGenomeSource) -I$(BASE_LIB_DIR)
	ar rcs lib$(LIB_NAME).a $(mapGenomeSource:%.cxx=%.o)


clean:
	rm *.o *.a
