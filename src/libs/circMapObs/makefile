cobsSource = circMapObs.cxx

BASE_LIB_DIR = ../

LIB_NAME = circMapObs


CXX = g++
CFLAGS = -c -g


cobs : $(cobsSource)
	$(CXX) $(CFLAGS) $(cobsSource) -I$(BASE_LIB_DIR)
	ar rcs lib$(LIB_NAME).a $(cobsSource:%.cxx=%.o)



clean:
	rm *.o *.a

