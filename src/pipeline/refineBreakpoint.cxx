/*

refineBreakpoint.cxx

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
#include "GFF/GFFtranscripts_utils.hxx"

using namespace std;


#define ARG_BASE_NAME 1
#define ARG_FILE_BLAT_POS 2
#define ARG_FILE_CS_POS 3
#define ARG_FILE_SAM 4

#define ARG_FILE_GENOME 5
#define ARG_FILE_GFF 6
#define ARG_FILE_GFF_TYPE 7
#define ARG_FILE_GFF_KEY 8



#define ARG_FILE_READS 9
#define ARG_FILE_MATES 10


int mod(int a, int b);

int main(int argc, char *argv[]){

  // Parse sam output and for each read with XM:i:1 try to refine the alignment and output the new version of the read + update the BLAT position file
  // output 2 sam files? 1 with un-modified sam output and 1 with refined
  // Modify BLAT OUTPUT positions


  char fname[8000];
  sprintf(fname, "%s-DEBREAKED.refined.fastq", argv[ARG_BASE_NAME]);
  FILE * nfDebreakedRefined = fopen(fname, "w");
  sprintf(fname, "%s-blat.map.refined", argv[ARG_BASE_NAME]);
  FILE * nfBlatRefined = fopen(fname, "w");
  sprintf(fname, "%s-XMzero.sam", argv[ARG_BASE_NAME]);
  FILE *nfSamXMzero = fopen(fname, "w");


  if(argc<2){
    fprintf(stderr, "%s: refines sam reads.\n", argv[0]);
    fprintf(stderr, "Usage:\n%s blat-mapping.pos lib-CS.pos map.sam ref-genome.fa exons.gff GFF key reads.fastq [mates.fastq]\n", argv[0]);
    return 1;
  }

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

    if( readSamRes < 0 ) // Skipping header lines
      continue;

    // Keep only uniquely matched reads
    int iFlag = sr.getFlag();
    if( iFlag != 0 && iFlag != 16 )
      continue;

    int nbMM = sr.getNbMM();
    if( nbMM > 1 )
      continue;

    // Locating the corresponding consensus
    string readID = sr.getQueryName();
    string chrID = sr.getRefName();

    //if( readID=="HSQ-7001360:223:H3L5LBCXX:1:2215:2124:95270")
    //fprintf(stderr, "BREAK!\n");

    CS_consensus cs;
    int csLg = -1;
    int alnLg = -1;
    int breakPos = -1;

    int res = naive_utils_getDebreakedCs(cs, readID, ficBlatPos, lineBlat, ficCsPos, ficReads, ficMates, paired, tab_cpt, &csLg, &alnLg, &breakPos);
    if( res != 0 ){
      fprintf(stderr, "ABORTING PROGRAM !\n");
      return -1;
    }

    if( sr.readReverseStrand() == true ){
      cs.reverseCpt(tab_cpt);
    }

    if( nbMM == 0 ){
      // Output samRead with no modification
      fprintf(nfSamXMzero, "%s\n", lineSam.c_str());
      //fprintf(nfBlatzero, "%s\n", lineBlat.c_str());
    } else {
      if(nbMM == 1){
	// Refine the alignment
	map<int,int> mc;
	if( sr.getGenomicCoords(mc) == false )
	  continue;

	vector<GFFtranscript> vt;
	getAllTranscriptsAtPosition(vt, mt, chrID, sr.getPos());

	int offset = 0;
	int res = cs.refine(&offset, sr, vt, mg);
	if( res == -1 ){
	  fprintf(stderr, "PROBLEM WITH %s\n", readID.c_str());
	}

	if( sr.readReverseStrand() == true ){
	  cs.reverseCpt(tab_cpt);
	  offset = 0 - offset;
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

      }
    }
  }

  fclose(nfSamXMzero);
  fclose(nfDebreakedRefined);
  fclose(nfBlatRefined);

  return 0;
}


int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}
