/*

  This program reads a map of observations and adds the strand information.

 */

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <vector>

#include "CircSeqConsensus/CS_consensus.hxx"
#include "UTILS_FASTA/utils_fasta.hxx"
#include "samRead/samRead.hxx"
#include "naive_utils.hxx"
#include "circMapCov-simple/circMapCovSimple.hxx"
#include "circMapIndels/circMapIndels.hxx"
#include "circMapObs/circMapObs.hxx"
#include "GFF/GFFtranscripts_utils.hxx"
#include "GFF/mapOfStrand.hxx"
#include "circMapGenome/circMapGenome.hxx"


using namespace std;

//#define ARG_FILE_DEBREAKED 1
#define ARG_FILE_MAP_OBS 1

#define ARG_FILE_GENOME 2
#define ARG_FILE_GFF 3
#define ARG_FILE_GFF_TYPE 4
#define ARG_FILE_GFF_KEY 5

int main(int argc, char *argv[]){

  if(argc<2){
    fprintf(stderr, "%s: adds strand information to map of observations.\n", argv[0]);
    fprintf(stderr, "Usage:\n%s mapObs.tab ref-genome.fa exons.gtf [GFF/GTF] gffKey\n", argv[0]);
    return 1;
  }

  int i;

  char tabBases[4] = {'A', 'T', 'C', 'G'};
  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);

  map<string, string> mg;
  fprintf(stderr, "Loading genome in memory (from '%s') ...", argv[ARG_FILE_GENOME]);
  loadFasta(mg, argv[ARG_FILE_GENOME]);
  fprintf(stderr, " done.\n");

  circMapGenome genomeMap;

  string gffType(argv[ARG_FILE_GFF_TYPE]);
  string gffKey(argv[ARG_FILE_GFF_KEY]);

  map<string, GFFtranscript> mt;
  vector<string> vSources;
  vector<string> vFeatures;
  vector<string> vChr;

  fprintf(stderr, "Loading annotations in memory (form %s)...", argv[ARG_FILE_GFF]);
  getAllTranscripts(argv[ARG_FILE_GFF], 0, mt, gffType, gffKey, 1, mg, tab_cpt, vSources, vFeatures, vChr);
  fprintf(stderr, " done.\n");

  fprintf(stderr, "Initializing the map of strands...");
  mapOfStrand mos(mg, mt);
  fprintf(stderr, " done.\n");

  string line, chrID, stmp;
  ifstream fic(argv[ARG_FILE_MAP_OBS]);
  getline(fic, line);
  while( getline(fic, line) ){
    istringstream iss(line);
    getline(iss, chrID, '\t');
    getline(iss, stmp, '\t');
    long gPos = atol(stmp.c_str());
    getline(iss, stmp, '\t');
    long nbObs = atol(stmp.c_str());
    char cStrand = mos.m_[chrID].at(gPos);
    fprintf(stdout, "%s\t%c\n", line.c_str(), cStrand);
  }


  /*
  fprintf(stderr, "Loading values into map of observations (from '%s') ...", argv[ARG_FILE_MAP_OBS]);
  genomeMap.loadObservations(argv[ARG_FILE_MAP_OBS], 1);
  fprintf(stderr, "done.\n");

  map<string, map<long,circMapGenomePos> >::iterator itc;
  for(itc = genomeMap.mg_.begin() ; itc!=genomeMap.mg_.end() ; itc++){
    string chrID = itc->first;
    map<long,circMapGenomePos>::iterator itp;
    for(itp=(itc->second).begin() ; itp!=(itc->second).end() ; itp++){
      long gPos = itp->first;
      long nbObs = (itp->second).getNbObservations();
      char cStrand = mos.m_[chrID].at(gPos);
      fprintf(stdout, "%s\t%d\t%d\t%c\n", chrID.c_str(), gPos, nbObs, cStrand);
    }
  }
  */

  return 0;
}
