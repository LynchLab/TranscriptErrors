/*

refieBreakpoint-withIndels.cxx

 */

#define _DEBOG_ 0
#undef _DEBOG_


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
#include "GFF/GFFtranscripts_utils.hxx"


using namespace std;


#define ARG_BASE_NAME 1
#define ARG_FILE_BLAT_POS 2
#define ARG_FILE_CS_POS 3

#define ARG_FILE_SAM 4
#define ARG_KEEP_ONLY_UNIQUE_MAPPING 5

#define ARG_FILE_GENOME 6

#define ARG_FILE_GFF 7
#define ARG_FILE_GFF_TYPE 8
#define ARG_FILE_GFF_KEY 9

#define ARG_FILE_READS 10
#define ARG_FILE_MATES 11


int mod(int a, int b);

int main(int argc, char *argv[]){

  // Parse sam output and for each read with XM:i:1 try to refine the alignment and output the new version of the read + update the BLAT position file
  // output 2 sam files? 1 with un-modified sam output and 1 with refined
  // Modify BLAT OUTPUT positions
  if(argc<2){
    fprintf(stderr, "%s: refines sam reads.\n", argv[0]);
    fprintf(stderr, "Usage:\n%s base_name blat-mapping.pos lib-CS.pos map.sam keepOnlyUniqueMapping [T/F] ref-genome.fa exons.gff gffType gffKey reads.fastq [mates.fastq]\n", argv[0]);
    return 1;
  }

  char fname[8000];
  sprintf(fname, "%s-DEBREAKED.refined.fastq", argv[ARG_BASE_NAME]);
  FILE * nfDebreakedRefined = fopen(fname, "w");
  sprintf(fname, "%s-blat.map.refined", argv[ARG_BASE_NAME]);
  FILE * nfBlatRefined = fopen(fname, "w");
  sprintf(fname, "%s-NMzero.sam", argv[ARG_BASE_NAME]);
  FILE *nfSamNMzero = fopen(fname, "w");

  sprintf(fname, "%s-LOST-CAUSE.sam", argv[ARG_BASE_NAME]);  // File for reads with edit size > 1
  FILE *nfLostCause = fopen(fname, "w");

  char tabBases[4] = {'A', 'T', 'C', 'G'};
  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);

  map<string, string> mg;
  fprintf(stderr, "Loading genome in memory (from '%s') ...", argv[ARG_FILE_GENOME]);
  loadFasta(mg, argv[ARG_FILE_GENOME]);
  fprintf(stderr, " done.\n");

  string gffType(argv[ARG_FILE_GFF_TYPE]);
  string gffKey(argv[ARG_FILE_GFF_KEY]);

  map<string, GFFtranscript> mt;
  vector<string> vSources;
  vector<string> vFeatures;
  vector<string> vSeqNames;
  getAllTranscripts(argv[ARG_FILE_GFF], 0, mt, gffType, gffKey, 1, mg, tab_cpt, vSources, vFeatures, vSeqNames);

  bool keepOnlyUniqueMapping = false;
  if( argv[ARG_KEEP_ONLY_UNIQUE_MAPPING][0]=='T' ){
    keepOnlyUniqueMapping = true;
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
  ifstream ficBlatPos(argv[ARG_FILE_BLAT_POS]);
  ifstream ficSam(argv[ARG_FILE_SAM]);

  string lineBlat, lineSam, stmp;
  samRead sr;

  int nbSamReadDone = 0;
  int readSamRes = 0;
  // Looping through the mapped reads (.sam file)
  while( readSamRes = sr.readNext(ficSam, lineSam) ){

#ifdef _DEBOG_
    string qName = sr.getQueryName();
    //    if(qName=="HSQ-7001360:381:HW2M5BCXY:2:1104:19185:43283" 
    //   || qName=="HSQ-7001360:381:HW2M5BCXY:2:1105:7226:66994"
    //  || qName=="HSQ-7001360:381:HW2M5BCXY:2:1101:9504:43981"
    // )
    //continue;
    fprintf(stderr, "%s\n", qName.c_str());
    fflush(stderr);
#endif

    if( readSamRes < 0 ) // Skipping header lines
      continue;

    // Keep only uniquely matched reads
    int iFlag = sr.getFlag();
    if( iFlag != 0 && iFlag != 16 )
      continue;

    if( keepOnlyUniqueMapping==true && sr.mOpt_["NH"]!="1" )
      continue;

    // Locating the corresponding consensus
    string readID = sr.getQueryName();
    string chrID = sr.getRefName();

    int eDist = sr.getEditDistance();
    if( eDist > 1 ){
      fprintf(nfLostCause, "%s\n", lineSam.c_str());
      continue;
    }

    if( eDist<0 ){
      fprintf(stderr, "WARNING: eDist<0 - readID:'%s'\n", readID.c_str());
    }

    CS_consensus cs;
    int csLg = -1;
    int alnLg = -1;
    int breakPos = -1;

    int res = naive_utils_getDebreakedCs(cs, readID, ficBlatPos, lineBlat, ficCsPos, ficReads, ficMates, paired, tab_cpt, &csLg, &alnLg, &breakPos);
    if( res < 0 ){
      fprintf(stderr, "ABORTING PROGRAM !\n");
      return -1;
    }
    if(res==1){
      fprintf(stderr, "WARNING, STRANGE CASE (%s)\n", readID.c_str());
      continue;
    }

    if( sr.readReverseStrand() == true ){
      cs.reverseCpt(tab_cpt);
    }

    if( eDist == 0 ){
      // Output samRead with no modification
      fprintf(nfSamNMzero, "%s\n", lineSam.c_str());
      //fprintf(nfBlatzero, "%s\n", lineBlat.c_str());
    } else {
      if(eDist == 1){
	// Refine the alignment

	vector<GFFtranscript> vt;
	getAllTranscriptsAtPosition(vt, mt, chrID, sr.getPos());

	int offset = 0;
	int res = cs.refine(&offset, sr, vt, mg);
	//fprintf(stderr, "RES from refine: %d (offset:%d)\n", res, offset);
	if( res == -1 ){
	  fprintf(stderr, "PROBLEM DURING REFINE OF %s\n", readID.c_str());
	}

	if( sr.readReverseStrand() == true ){
	  cs.reverseCpt(tab_cpt);
	  offset = 0 - offset; // DEBOG!!!
	}

	int newBreakPos = mod( (breakPos + offset), csLg); // CHECK HERE THAT I DID NOT FORGET ANY SPECIFIC CASE.

	// Output the new -DEBREAKED.fastq file
	string refinedCs = cs.getConsensus();
	string refinedQual = cs.getCsQuality('"', 'J');
	fprintf(nfDebreakedRefined, 
		"@%s\n%s\n%c\n%s\n",
		readID.c_str(),
		refinedCs.c_str(),
		'+',
		refinedQual.c_str()
		);

	//cs.prettyPrint(nfDebreakedRefined, true);

	// Update the blat-mapping.size file
	fprintf(nfBlatRefined, "%s\t%d\t%d\t%d\n", readID.c_str(), csLg, alnLg, newBreakPos);
	//fflush(nfBlatRefined);
	//fprintf(stderr, "WRITING: %s\t%d\t%d\t%d\n", readID.c_str(), csLg, alnLg, newBreakPos);
	//fflush(stderr);
      }
    }
  }

  fclose(nfSamNMzero);
  fclose(nfDebreakedRefined);
  fclose(nfBlatRefined);

  return 0;
}


int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}
