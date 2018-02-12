/*

testMap-indels.cxx

 */

#include "circMapIndels/circMapIndels.hxx"

#define ARG_FILE_GENOME 1
#define ARG_FILE_SAM_LIST 2


using namespace std;


int main(int argc, char *argv[]){

  if(argc<2){
    fprintf(stderr, "%s: builds a map of indel positions (with coverage) from a list of sam files.\n", argv[0]);
    fprintf(stderr, "Usage:\n%s genome.fa listSam.txt [text file with name of one sam file per line]\n", argv[0]);
    return 1;
  }

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

  fprintf(stderr, "Initializing map of indels ...");
  fflush(stderr);
  C_circMapIndels mapi(argv[ARG_FILE_GENOME]);
  fprintf(stderr, "done.\n");
  fflush(stderr);

  //mapc.printSize(stdout);
  //fflush(stdout);

  string line;
  ifstream fic(argv[ARG_FILE_SAM_LIST]);
  while( getline(fic, line) ){
    fprintf(stderr, "%s\n", line.c_str());
    fflush(stderr);
    mapi.update((char *)line.c_str());
  }


  mapi.print(stdout);

  return 0;
}
