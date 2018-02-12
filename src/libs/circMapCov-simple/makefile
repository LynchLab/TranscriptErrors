mapSource = circMapCovSimple.cxx

BASE_LIB_DIR = ../

LIB_NAME = circMapCovSimple

CXX = g++
CFLAGS = -c -g

all: map

map: $(mapSource)
	$(CXX) $(CFLAGS) $(mapSource) -I$(BASE_LIB_DIR)
	ar rcs lib$(LIB_NAME).a $(mapSource:%.cxx=%.o)


clean:
	rm *.o *.a


