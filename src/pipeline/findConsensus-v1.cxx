/**

@file
@author
@date
@brief This program finds consensus (=repeats) in RNAseq reads generated with the circ-seq method.

Finding repeats is based on a very simple algorithm: searching for the second occurence of the first N (typically 30) nucleotides in the sequence (allowing for some mismatchs).


*/

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include "CircSeqConsensus/CS_consensus.hxx"
#include "UTILS_FASTA/utils_fasta.hxx"
#include "naive_utils.hxx"


#define ARG_FILE_BASENAME 1
#define ARG_FILE_READS 2
#define ARG_FILE_MATES 3


using namespace std;



int main(int argc, char *argv[]){

  if( argc < 2 ){
    fprintf(stderr, "%s: finds consensus sequences from circseq-generated data.\n", argv[0]);
    fprintf(stderr, "Usage:\n%s base_name <reads1.fastq,reads2.fastq,...> <mates1.fastq,mates2.fastq,...>\n", argv[0]);
    fprintf(stderr, "[base_name] is used to write 2 output files: base_name-consensus.fa and base_name-NoConsensus.fa representing the reads with and without a consensus respectively.\n");
    return 1;
  }

  char CHAR_SEP = ',';
  int minQual = '"';
  int maxQual = 'J';

  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);

  int minRepeatSize = 30;
  int maxMM = 2;
  double minID = 0.9;

  string line, readID, readSeq, readQual, mateID, mateSeq, mateQual, stmp;
  int nbHaveConsensus = 0;
  int nbTotal = 0;

  bool paired = false;
  vector<string> vReads;
  vector<string> vMates;

  string baseName(argv[1]);
  string sReads(argv[2]);

  fillVectorFromArgs(vReads, sReads, CHAR_SEP); // Loading the list of file names for (left)-reads


  //***************************************************************/
  // Preparing the output files
  char ctmp[8000];
  sprintf(ctmp, "%s-DCS.fastq", argv[ARG_FILE_BASENAME]);
  FILE *nfCS = fopen(ctmp, "w");

  sprintf(ctmp, "%s-CS.pos", argv[ARG_FILE_BASENAME]);
  FILE *nfPos = fopen(ctmp, "w");

  FILE *nfNCS, *nfNCS_left, *nfNCS_right;
  if( argc>ARG_FILE_MATES ){
    paired = true;
    string sMates(argv[3]);
    fillVectorFromArgs(vMates, sMates, CHAR_SEP); // Loading the list of file names for right-reads (if paired-end sequencing)
    sprintf(ctmp, "%s-NCS-R1.fastq", argv[ARG_FILE_BASENAME]);
    nfNCS_left = fopen(ctmp, "w");
    sprintf(ctmp, "%s-NCS-R2.fastq", argv[ARG_FILE_BASENAME]);
    nfNCS_right = fopen(ctmp, "w");
  } else {
    sprintf(ctmp, "%s-NCS.fastq", argv[ARG_FILE_BASENAME]);
    nfNCS = fopen(ctmp, "w");
  }
  // Output files are ready to be used
  ////////////////////////////////////////////////////////////////////

  ifstream ficReads, ficMates;
  fastqRead read;
  fastqRead mate;

  int nbFiles = vReads.size();
  for(int iFile=0 ; iFile<nbFiles ; iFile++){
    ficReads.open(vReads[iFile].c_str(), std::ifstream::in);
    if( paired==true ){
      ficMates.open(vMates[iFile].c_str(), std::ifstream::in);
    }

    while( read.readNext(ficReads) ){ // Reading the reads sequences one by one from the fastq file
      if( paired==true && !mate.readNext(ficMates) ){
	fprintf(stderr, "ERROR: '%s' file ended before '%s'\nTERMINATING PROGRAM.\n", vReads[iFile].c_str(), vMates[iFile].c_str());
	exit(2);
      }

      if( paired==true && read.id_ != mate.id_ ){
	fprintf(
		stderr, 
		"ERROR: '%s' and '%s' are de-synchronized (readID:'%s' - mateID:'%s')\nTERMINATING PROGRAM.\n", 
		vReads[iFile].c_str(), vMates[iFile].c_str(), read.id_.c_str(), mate.id_.c_str()
		);
      }

      read.trimLowQuality('#');
      CS_consensus cs;
      if( paired==true ){ 
	mate.trimLowQuality('#');
	string mateSeq = my_reverse_cpt(mate.seq_, tab_cpt);
	string mateQual = my_reverse(mate.squal_);
	cs.make(read.seq_, read.squal_, mateSeq, mateQual, minRepeatSize, maxMM, minID);
      } else {
	cs.make(read.seq_, read.squal_, minRepeatSize, maxMM, minID);
      }

      string scs = cs.getConsensus();
      if(scs.size() >= minRepeatSize ){

	fprintf(stdout, ">%s\n", read.id_.c_str());
	cs.prettyPrint(stdout, true);

	string sQual = cs.getCsQuality(minQual, maxQual);
	fprintf(nfCS, "@%s\n%s%s\n%c\n%s%s\n", read.id_.c_str(), scs.c_str(), scs.c_str(), '+', sQual.c_str(), sQual.c_str());
	fprintf(nfPos, "%s\t%d", read.id_.c_str(), cs.getPosStartInRead());
	if( paired==true ){
	  fprintf(nfPos, "\t%d", cs.getPosStartInMate());
	}
	fprintf(nfPos, "\n");
	nbHaveConsensus++;
      } else {
	if( paired==true ){ 
	  read.fullPrint(nfNCS_left);
	  mate.fullPrint(nfNCS_right);
	} else {
	  read.fullPrint(nfNCS);
	}
      }

    nbTotal++;

    }//END of loop all entries in a give (pair of) file(s)

    ficReads.close();
    if( paired==true ){
      ficMates.close();
    }

  }//END of loop all files


  fclose(nfCS);
  fclose(nfPos);
  if( paired==true ){
    fclose(nfNCS_left);
    fclose(nfNCS_right);
  } else {
    fclose(nfNCS);
  }

  cerr<<nbHaveConsensus
      <<" out of "
      <<nbTotal
      <<" reads have a consensus."
      <<endl;

  return 0;
}



