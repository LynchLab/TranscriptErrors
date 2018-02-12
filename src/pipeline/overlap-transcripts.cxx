#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <vector>

#include "UTILS_FASTA/utils_fasta.hxx"
#include "GFF/GFFtranscripts_utils.hxx"


#define ARG_FILE_GENOME 1
#define ARG_FILE_GFF 2


using namespace std;


int main(int argc, char *argv[]){



  char tabBases[4] = {'A', 'T', 'C', 'G'};
  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);


  //map<string, string> mg;
  //fprintf(stderr, "Loading genome in memory (from '%s') ...", argv[ARG_FILE_GENOME]);
  //loadFasta(mg, argv[ARG_FILE_GENOME]);
  //fprintf(stderr, " done.\n");


  map<string, GFFtranscript> mt;
  vector<string> vSources;
  vector<string> vFeatures;
  getAllTranscripts(argv[ARG_FILE_GFF], 0, mt, "GTF", "transcript_id", 1, vSources, vFeatures);
 
  //fprintf(stderr, "map size: %d\n", (int)mt.size());

  map<string, GFFtranscript>::iterator it;
  for(it=mt.begin() ; it!=mt.end() ; it++){
    GFFtranscript gt = it->second;
    long posMin, posMax;
    gt.getGenomicMinMaxPos(&posMin, &posMax);
    fprintf(stdout, "%s\t%ld\t%ld\n", gt.id_.c_str(), posMin, posMax);
  }


  return 0;
}
