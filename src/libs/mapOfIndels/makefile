mapSource = indelCoverage.cxx mapOfIndel.cxx

BASE_LIB_DIR = ../

LIB_NAME = mapOfIndel


CXX = g++
CFLAGS = -c -g


map : $(mapSource)
	$(CXX) $(CFLAGS) $(mapSource) -I$(BASE_LIB_DIR)
	ar rcs lib$(LIB_NAME).a $(mapSource:%.cxx=%.o)


clean:
	rm *.o *.a

