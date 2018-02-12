/**

@file
@author
@date
@brief Builds the map of observed 'genotypes' from a list of mapped reads (sam files)

This program takes two files as input: the genome in fasta format and a file containing a list of sam files (one sam file per line).
For each position in the genome, the program outputs the number of A, T, C, and G observations.
This is quite similar to what mipleup would do and this step in the pipeline could actually be substituted by using mpileup.

Output: a file containing one line per position in the genome with the number of A, T, C, and G observations.

 */

#include "circMapCov-simple/circMapCovSimple.hxx"

#define ARG_FILE_GENOME 1
#define ARG_FILE_SAM_LIST 2


using namespace std;


int main(int argc, char *argv[]){

  int i;

  char tabBases[4] = {'A', 'T', 'C', 'G'};
  int tabCor['z' + 'Z'];
  for(i=0 ; i<'z'+'Z' ; i++){
    tabCor[i] = -1;
  }

  for(i=0 ; i<4 ; i++){
    tabCor[tabBases[i]] = i;
    tabCor[tolower(tabBases[i])] = i;
  }

  fprintf(stderr, "Initializing map ...");
  fflush(stderr);
  circMapCovSimple mapc(argv[ARG_FILE_GENOME]);
  fprintf(stderr, "done.\n");
  fflush(stderr);

  //mapc.printSize(stdout);
  //fflush(stdout);

  string line;
  ifstream fic(argv[ARG_FILE_SAM_LIST]);
  while( getline(fic, line) ){
    fprintf(stderr, "%s\n", line.c_str());
    fflush(stderr);
    mapc.update((char *)line.c_str(), tabCor, 30);
  }


  mapc.print(stdout, tabBases);

  return 0;
}
