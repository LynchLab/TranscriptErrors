#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <vector>

#include "UTILS_FASTA/utils_fasta.hxx"
#include "GFF/GFFtranscripts_utils.hxx"
#include "circMapObs/circMapObs.hxx"


#define ARG_FILE_FILESOBS 1
#define ARG_FILE_GENOME 2
#define ARG_FILE_GFF 3


using namespace std;


int main(int argc, char *argv[]){

  string line;

  char tabBases[4] = {'A', 'T', 'C', 'G'};
  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);


  map<string, string> mg;
  fprintf(stderr, "Loading genome in memory (from '%s') ...", argv[ARG_FILE_GENOME]);
  loadFasta(mg, argv[ARG_FILE_GENOME]);
  fprintf(stderr, " done.\n");


  map<string, GFFtranscript> mt;
  vector<string> vSources;
  vector<string> vFeatures;
  getAllTranscripts(argv[ARG_FILE_GFF], 0, mt, "GTF", "transcript_id", 1, mg, tab_cpt, vSources, vFeatures);
 
  //fprintf(stderr, "map size: %d\n", (int)mt.size());

  vector<string> vFiles;
  ifstream ficFilesObs(argv[ARG_FILE_FILESOBS]);
  while( getline(ficFilesObs, line) ){
    vFiles.push_back(line);
  }
  circMapObs mobs(vFiles);

mobs.print(stdout);

  map<string, GFFtranscript>::iterator it;
  for(it=mt.begin() ; it!=mt.end() ; it++){
    GFFtranscript gt = it->second;
    long posMin, posMax;
    gt.getGenomicMinMaxPos(&posMin, &posMax);
//fprintf(stdout, ">%s\t%s\n", gt.id_.c_str(), gt.seq_.c_str());
//fprintf(stdout, "%s\n", gt.seq_.c_str());
  }


  return 0;
}
