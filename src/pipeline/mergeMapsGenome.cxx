/**

@file
@author
@date
@brief Merges multiple maps of observed 'genotypes'

This program takes one file as input: a text file containing the path to  map files, one per line.

Output: a file containing one line per position in the genome with the number of A, T, C, and G observations summed across all the files in the input.

 */

#include <stdlib.h>

#include "circMapGenome/circMapGenome.hxx"
#include "circMapCov-simple/circMapCovSimple.hxx"


//#define ARG_FILE_GENOME 1
#define ARG_FILE_LIST 1
#define ARG_FILE_NUCLEOTIDES_OUT 2
#define ARG_FILE_INDELS_OUT 3


using namespace std;


int main(int argc, char *argv[]){

  if(argc<3){
    fprintf(stderr, "%s: merges multiple maps of coverage and indels from a list of map files.\n", argv[0]);
    fprintf(stderr, "Usage:\n%s genome.fa files.list \n", argv[0]);
    return 1;
  }

  int i;


  FILE * nf_mapCoverage = fopen(argv[ARG_FILE_NUCLEOTIDES_OUT], "w");
  FILE * nf_mapIndels;
  if( argc>ARG_FILE_NUCLEOTIDES_OUT ){
   nf_mapIndels = fopen(argv[ARG_FILE_INDELS_OUT], "w");
  }


  char tabBases[4] = {'A', 'T', 'C', 'G'};
  int tabCor['z' + 'Z'];
  for(i=0 ; i<'z'+'Z' ; i++){
    tabCor[i] = -1;
  }

  for(i=0 ; i<4 ; i++){
    tabCor[tabBases[i]] = i;
    tabCor[tolower(tabBases[i])] = i;
  }


  //map<string,string> mg;
  //loadFasta(mg, argv[ARG_FILE_GENOME]);

  fprintf(stderr, "Initializing map ...");
  fflush(stderr);
  //circMapGenome mapGenome(mg);
  circMapGenome mapGenome;
  fprintf(stderr, "done.\n");
  fflush(stderr);

  string line, sname;
  bool hasIndels = false;

  ifstream fic(argv[ARG_FILE_LIST]);
  while( getline(fic, line) ){
    istringstream iss(line);
    getline(iss, sname, '\t');
    mapGenome.loadMapCoverage((char *)sname.c_str(), true);
    sname = "";
    getline(iss, sname, '\t');
    if( sname.size() > 1 ){
      mapGenome.loadMapIndels((char *)sname.c_str(), true);
      hasIndels = true;
    }
  }

  mapGenome.printMapCoverage(nf_mapCoverage, tabBases);
  if( hasIndels==true ){
    mapGenome.printMapIndels(nf_mapIndels);
  }

  return 0;
}
