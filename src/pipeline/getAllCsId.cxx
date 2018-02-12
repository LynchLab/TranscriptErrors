/*

  Just a utility program to get the distribution of consensus prct identity.

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



using namespace std;

//#define ARG_FILE_DEBREAKED 1
#define ARG_FILE_BLAT_POS 1
#define ARG_FILE_CS_POS 2
#define ARG_FILE_SAM 3
#define ARG_KEEP_ONLY_UNIQUE 4

#define ARG_FILE_READS 5
#define ARG_FILE_MATES 6



int main(int argc, char *argv[]){

  if(argc<2){
    fprintf(stderr, "%s: computes the prct identity for every consensus.\n", argv[0]);
    fprintf(stderr, "Usage:\n%s blat-mapping.pos lib-CS.pos map.sam keepOnlyUniqueMapped [T/F] reads.fastq [mates.fastq]\n", argv[0]);
    return 1;
  }

  int i;

  char tabBases[4] = {'A', 'T', 'C', 'G'};
  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);

  bool keepOnlyUniqueMapped = false;
  if( argv[ARG_KEEP_ONLY_UNIQUE][0]=='T' ){
    keepOnlyUniqueMapped = true;
  }

  bool paired = false;
  fprintf(stderr, "PAIRED FILES?");
  if( argc>ARG_FILE_MATES ){
    paired = true;
    fprintf(stderr, " YES\n");
  } else {
    fprintf(stderr, " NO.\n");
  }

  ifstream ficReads, ficMates;
  ficReads.open(argv[ARG_FILE_READS]);
  if( paired==true ){
    ficMates.open(argv[ARG_FILE_MATES]);
  }

  ifstream ficCsPos(argv[ARG_FILE_CS_POS]);

  string lineBlat, lineSam, stmp;

  ifstream ficBlatPos(argv[ARG_FILE_BLAT_POS]);
  ifstream ficSam(argv[ARG_FILE_SAM]);

  samRead sr;
  CS_consensus cs;


  int nbSamReadDone = 0;
  int resFromSam = 0;
  // Looping through the mapped reads (.sam file)
  while( resFromSam=sr.readNext(ficSam) ){

    if(resFromSam < 0) // skipping header lines
      continue;

    // Keep only uniquely matched reads
    int iFlag = sr.getFlag();
    if( iFlag != 0 && iFlag != 16 )
      continue;

    if( keepOnlyUniqueMapped==true && sr.mOpt_["NH"]!="1" )
      continue;

    ///////////////////////////////////////////////////////////////////
    //                DEBOG ONLY - REMOVE FOR PRODUCTION 
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( sr.isOnIntron()==true )
      continue;

    // Locating the corresponding consensus
    string queryName = sr.getQueryName();
    cs.clean();
    // Building the consensus
    int res = naive_utils_getDebreakedCs(cs, queryName, ficBlatPos, ficCsPos, ficReads, ficMates, paired, tab_cpt);

    if( res<0 ){
      fprintf(stderr, "ABORTING PROGRAM.\n");
      return 1;
    }
    
    if( res>0 )
      continue;
    
    // Here, look at CS
    if( sr.readReverseStrand() == true ){
      cs.reverseCpt(tab_cpt);
    }
    
    fprintf(stdout, "%s\t%.6f\n", queryName.c_str(), cs.getPrctId());
    nbSamReadDone++;

    if(nbSamReadDone%100000 == 0){
      fprintf(stderr, "%d sam reads done.\n", nbSamReadDone);
    }
  }
  
  return 0;
}



