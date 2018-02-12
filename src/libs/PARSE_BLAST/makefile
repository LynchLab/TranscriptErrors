PB = parse_blast

CXX = g++
CFLAGS = -c

all: pb

pb: $(PB).cxx
	$(CXX) $(CFLAGS) $(PB).cxx -o $(PB).o
	ar rcs lib$(PB).a $(PB).o

clean:
	rm *.o *.a
