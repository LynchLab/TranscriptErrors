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

#include <stdlib.h>

#include "circMapGenome/circMapGenome.hxx"
#include "circMapCov-simple/circMapCovSimple.hxx"


//#define ARG_FILE_GENOME 1
#define ARG_FILE_SAM_LIST 1
#define ARG_MIN_SCORE 2

#define ARG_FILE_GENOME 3
//#define ARG_FILE_GFF 4
//#define ARG_FILE_GFF_TYPE 5
//#define ARG_FILE_GFF_KEY 6

#define ARG_FILE_MAP_COVERAGE_OUT 4
#define ARG_FILE_MAP_INDELS_OUT 5


using namespace std;


int main(int argc, char *argv[]){

  if(argc<3){
    fprintf(stderr, "%s: computes the maps of coverage and indels from a list of sam files.\n", argv[0]);
    //fprintf(stderr, "Usage:\n%s samFiles.list minScore genome.fa exons.gff gffType gffKey map.tab map-indels.tab\n", argv[0]);
    fprintf(stderr, "Usage:\n%s samFiles.list minScore genome.fa map.tab map-indels.tab\n", argv[0]);
    return 1;
  }

  int i;

  int minScore = atoi(argv[ARG_MIN_SCORE]);

  FILE * nf_mapCoverage = fopen(argv[ARG_FILE_MAP_COVERAGE_OUT], "w");
  FILE * nf_mapIndels = fopen(argv[ARG_FILE_MAP_INDELS_OUT], "w");


  char tabBases[4] = {'A', 'T', 'C', 'G'};
  int tabCor['z' + 'Z'];
  for(i=0 ; i<'z'+'Z' ; i++){
    tabCor[i] = -1;
  }

  for(i=0 ; i<4 ; i++){
    tabCor[tabBases[i]] = i;
    tabCor[tolower(tabBases[i])] = i;
  }

  /*
  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);

  fprintf(stderr, "Preparing a map of strands ...");
  string gffType(argv[ARG_FILE_GFF_TYPE]);
  string gffKey(argv[ARG_FILE_GFF_KEY]);
  map<string, GFFtranscript> mt;
  vector<string> vSources;
  vector<string> vFeatures;
  vector<string> vChr;
  getAllTranscripts(argv[ARG_FILE_GFF], 0, mt, gffType, gffKey, 1, mg, tab_cpt, vSources, vFeatures, vChr);
  mapOfStrand mos(mg, mt);
  fprintf(stderr, " done.\n");
  */

  fprintf(stderr, "Initializing map ...");
  fflush(stderr);
  circMapGenome mapGenome;
  fprintf(stderr, "done.\n");
  fflush(stderr);

  //mapc.printSize(stdout);
  //fflush(stdout);

  string line;
  samRead sr;

  long nbReadsDone = 0;

  ifstream fic(argv[ARG_FILE_SAM_LIST]);
  while( getline(fic, line) ){
    fprintf(stderr, "%s\n", line.c_str());
    ifstream ficSam(line.c_str());
    while( sr.readNext(ficSam) ){
      string chrID = sr.getRefName();
      map<string, map<long,circMapGenomePos> >::iterator itc = mapGenome.mg_.find(chrID);
      if( itc==mapGenome.mg_.end() ){
	map<long,circMapGenomePos> mtmp;
	itc = (mapGenome.mg_.insert( map<string, map<long,circMapGenomePos> >::value_type(chrID, mtmp) )).first;
      }

      //string qName = sr.getQueryName();
      //if(qName=="HSQ-7001360:296:HTCLJBCXX:2:2209:2673:82912"){
      //fprintf(stderr, "DEBOG.\n");
      //}

      mapGenome.processSamRead_forCoverage(sr, itc, tabCor, minScore);
      mapGenome.updateMapIndels(sr, itc);

      nbReadsDone++;
      if(nbReadsDone%100000 == 0){
	fprintf(stderr, "%d reads processed...\n", nbReadsDone);
      }

    }
  }

  map<string,string> mg;
  fprintf(stderr, "Loading genome in memory (from '%s') ...", argv[ARG_FILE_GENOME]);
  loadFasta(mg, argv[ARG_FILE_GENOME]);
  fprintf(stderr, " done.\n");

  fprintf(stderr, "Updating reference base ...");
  mapGenome.updateAllRefBase(mg);
  fprintf(stderr, " done.\n");

  fprintf(stderr, "Outputintg the map of coverage to: '%s' ...", argv[ARG_FILE_MAP_COVERAGE_OUT]);
  mapGenome.printMapCoverage(nf_mapCoverage, tabBases);
  fprintf(stderr, " done.\n");

  fprintf(stderr, "Outputintg the map of indels to: '%s' ...", argv[ARG_FILE_MAP_INDELS_OUT]);
  mapGenome.printMapIndels(nf_mapIndels);
  fprintf(stderr, " done.\n");

  fclose(nf_mapCoverage);
  fclose(nf_mapIndels);

  return 0;
}
