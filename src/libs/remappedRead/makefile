remappedRead_source = remappedRead.cxx remappedRead_utils.cxx

BASE_LIB_DIR = ../

LIB_NAME = remappedRead

CXX = g++
CFLAGS = -g -c


remappedRead : $(remappedRead_source)
	$(CXX) $(CFLAGS) $(remappedRead_source) -I$(BASE_LIB_DIR)
	ar rcs lib$(LIB_NAME).a $(remappedRead_source:%.cxx=%.o)


clean:
	rm *.o *.a
