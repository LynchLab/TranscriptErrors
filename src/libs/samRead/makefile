BASE_LIB_DIR = ../

SAMREAD = samRead

CXX = g++
CFLAGS = -g -c



samRead : $(SAMREAD).cxx
	$(CXX) $(CFLAGS) $(SAMREAD).cxx -I$(BASE_LIB_DIR) -o $(SAMREAD).o
	ar rcs lib$(SAMREAD).a $(SAMREAD).o

clean:
	rm *.o *.a
