/*

  This program adds the base to a map of observations.
  If the strand information is ambiguous, the base is '?'

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
#include "circMapObs/circMapObs.hxx"
#include "GFF/GFFtranscripts_utils.hxx"
#include "GFF/mapOfStrand.hxx"
#include "circMapGenome/circMapGenome.hxx"


using namespace std;


#define ARG_FILE_MAP 1
#define ARG_FILE_GENOME 2

#define ARG_FILE_GFF 3
#define ARG_GFF_TYPE 4
#define ARG_GFF_KEY 5


int main(int argc, char *argv[]){

  if(argc<4){
    fprintf(stderr, "%s: adds the base nucleotide information to a map of observations.\n", argv[0]);
    fprintf(stderr, "Usage:\n%s observations.tab genome.fa annotations.gtf [GFF/GTF] key [Parent, transcript_id, ...]\n", argv[0]);
    return 1;
  }

  string gffType(argv[ARG_GFF_TYPE]);
  string sKey(argv[ARG_GFF_KEY]);

  char tabBases[4] = {'A', 'T', 'C', 'G'};
  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);

  map<string, string> mg;
  fprintf(stderr, "Loading genome in memory (from '%s') ...", argv[ARG_FILE_GENOME]);
  loadFasta(mg, argv[ARG_FILE_GENOME]);
  fprintf(stderr, " done.\n");

  circMapGenome genomeMap;
  fprintf(stderr, "Initializating the map of observations ...");
  genomeMap.loadObservations(argv[ARG_FILE_MAP]);
  fprintf(stderr, " done.\n");

  map<string, GFFtranscript> mt;
  //vector<string> vSources;
  //vector<string> vFeatures;
  //vSources.push_back("RefSeq");
  //vSources.push_back("ena");
  //vSources.push_back("ensembl");

  //vFeatures.push_back("transcript");
  //vFeatures.push_back("exon");
  //vFeatures.push_back("gene");
  //vFeatures.push_back("rRNA");
  //vFeatures.push_back("tRNA");
  //vFeatures.push_back("ncRNA");

  getAllTranscripts(argv[ARG_FILE_GFF], 0, mt, gffType, sKey, 1);
  fprintf(stderr, "map size: %d\n", (int)mt.size());

  //getAllTranscripts(argv[ARG_FILE_GFF], 0, mt, "GFF", "Parent", 1);
  //getAllTranscripts(argv[ARG_FILE_GFF], 5, mt, "GTF", "gene_id", 1, mg, tab_cpt, vSources, vFeatures);
  //getAllTranscripts(argv[ARG_FILE_GFF], 5, mt, "GFF", "ID", 5, mg, tab_cpt, vSources, vFeatures);
  mapOfStrand mos(mg, mt);

  genomeMap.printObservationsWithBase(stdout, mg, mos, false);

  return 0;
}

