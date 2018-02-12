# Makefile for GFF library
BASE_LIB_DIR = ../

GFF = GFF

GFF_SOURCE = GFFentry.cxx GFFtranscripts.cxx GFFtranscripts_utils.cxx GFFutils.cxx mapOfStrand.cxx

CXX = g++
CFLAGS = -g -c


all: gff

gff: $(GFF_SOURCE)
	$(CXX) $(CFLAGS) $(GFF_SOURCE) -I$(BASE_LIB_DIR)
	ar rcs lib$(GFF).a $(GFF_SOURCE:%.cxx=%.o)

clean:
	rm *.o *.a
