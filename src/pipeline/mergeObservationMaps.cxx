/*

  Merges multiple files of observations (each input file contains 3 columns: chromosome position #observation / the output file 
  contains the sum of the third columns in all input files).

 */

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <vector>

#include "circMapObs/circMapObs.hxx"
#include "circMapGenome/circMapGenome.hxx"
//#include "GFF/mapOfStrand.hxx"
//#include "UTILS_FASTA/utils_fasta.hxx"

using namespace std;


#define ARG_FILE_FILES_OBS 1
#define ARG_OUTPUT_ZEROS 2

int main(int argc, char *argv[]){

  if(argc<2){
    fprintf(stderr, "%s: merges observation files.\n", argv[0]);
    fprintf(stderr, "Usage:\n%s files.list > merged.tab (files.list contains the list of files to merge, one per line)\n", argv[0]);
    return 1;
  }

  bool outputZeros = true;
  if( argc>ARG_OUTPUT_ZEROS && argv[ARG_OUTPUT_ZEROS][0]=='F' ){
    outputZeros = false;
  }

  circMapGenome genomeMap;

  string line;
  ifstream fic(argv[ARG_FILE_FILES_OBS]);
  while( getline(fic, line) ){
    fprintf(stderr, "Adding observations from '%s' ...", line.c_str());
    genomeMap.addObservations((char *)line.c_str(), 0);
    fprintf(stderr, " done.\n");
  }

  fprintf(stderr, "Outputing themerged map of observations ...");
  genomeMap.printObservations(stdout, outputZeros);
  fprintf(stderr, " done.\n");


  return 0;
}

