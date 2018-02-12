UF = utils_fasta

CXX = g++
CFLAGS = -c


all: uf


uf: $(UF).cxx
	$(CXX) $(CFLAGS) $(UF).cxx -o $(UF).o
	ar rcs lib$(UF).a $(UF).o

clean:
	rm *.o *.a
