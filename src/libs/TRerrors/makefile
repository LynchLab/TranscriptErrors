BASE_LIB_DIR = ../

TRE_SOURCE = TRerrors.cxx

TRE_LIB = TRerrors


CXX = g++
CFLAGS = -g -c


all: tre

tre: $(TRE_SOURCE)
	$(CXX) $(CFLAGS) $(TRE_SOURCE) -I$(BASE_LIB_DIR)
	ar rcs lib$(TRE_LIB).a $(TRE_SOURCE:%.cxx=%.o)

clean:
	rm *.o *.a
