/*

  The program that parses the sam output from the mapping of consensus sequences and calls candidate transcription errors.

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


using namespace std;


int processMappedCs(samRead &sr, CS_consensus &cs, circMapCovSimple &mapc, map<string, string> &mg, char *tabBases, char *tab_cpt, long *nbErrors, long *nbSearch, long *tabObs, mapOfStrand &mos, circMapObs &mobs, int minCov, double macFracAlter);
void printDetailedMM(FILE *nf, CS_consensus &cs, string &genomeSeq, int posCs);


//#define ARG_FILE_DEBREAKED 1
#define ARG_FILE_BLAT_POS 1
#define ARG_FILE_CS_POS 2
#define ARG_FILE_SAM 3
#define ARG_FILE_GENOME 4
#define ARG_FILE_GFF 5
#define ARG_FILE_GFF_TYPE 6
#define ARG_FILE_GFF_KEY 7

#define ARG_FILE_MAP 8

#define ARG_MIN_COV 9
#define ARG_MAX_FRAC_ALTER 10
#define ARG_FILE_OBS 11

#define ARG_FILE_READS 12
#define ARG_FILE_MATES 13



int main(int argc, char *argv[]){

  if(argc<2){
    fprintf(stderr, "%s: calls transcription errors from circseq data.\n", argv[0]);
    fprintf(stderr, "Usage:\n%s blat-mapping.pos lib-CS.pos map.sam ref-genome.fa ref.gtf [GFF/GTF] gffKey map.tab minCov maxFracAlter file-obs reads.fastq [mates.fastq]\n", argv[0]);
    return 1;
  }

  int i;

  long nbErrors = 0;
  long nbSearch = 0;

  int minCov = atoi(argv[ARG_MIN_COV]);
  double maxFracAlter = atof(argv[ARG_MAX_FRAC_ALTER]);

  char tabBases[4] = {'A', 'T', 'C', 'G'};
  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);

  long tabObs['z'+'Z'];
  for(i=0 ; i<'z'+'Z' ; i++){
    tabObs[i] = 0;
  }

  map<string, string> mg;
  fprintf(stderr, "Loading genome in memory (from '%s') ...", argv[ARG_FILE_GENOME]);
  loadFasta(mg, argv[ARG_FILE_GENOME]);
  fprintf(stderr, " done.\n");

  fprintf(stderr, "Initializating the map of observations ...");
  circMapObs mobs(mg);
  fprintf(stderr, " done.\n");

  string gffType(argv[ARG_FILE_GFF_TYPE]);
  string gffKey(argv[ARG_FILE_GFF_KEY]);

  map<string, GFFtranscript> mt;
  vector<string> vSources;
  vector<string> vFeatures;
  //vSources.push_back("RefSeq");
  //vSources.push_back("ena");
  //vSources.push_back("ensembl");


  //vFeatures.push_back("transcript");
  //vFeatures.push_back("exon");
  //vFeatures.push_back("gene");
  //vFeatures.push_back("rRNA");
  //vFeatures.push_back("tRNA");
  //vFeatures.push_back("ncRNA");

  //getAllTranscripts(argv[ARG_FILE_GFF], 5, mt, "GTF", "transcript_id", 1, mg, tab_cpt, vSources, vFeatures);
  //getAllTranscripts(argv[ARG_FILE_GFF], 5, mt, "GFF", "Parent", 1, mg, tab_cpt, vSources, vFeatures);
  //getAllTranscripts(argv[ARG_FILE_GFF], 5, mt, "GFF", "ID", 5, mg, tab_cpt, vSources, vFeatures);
  getAllTranscripts(argv[ARG_FILE_GFF], 0, mt, gffType, gffKey, 1, mg, tab_cpt, vSources, vFeatures);
  mapOfStrand mos(mg, mt);

  fprintf(stderr, "Initializating map of coverage ...");
  circMapCovSimple mapc(mg);
  fprintf(stderr, " done.\n");
  fprintf(stderr, "Loading values into map of coverage (from '%s') ...", argv[ARG_FILE_MAP]);
  mapc.load(argv[ARG_FILE_MAP]);
  fprintf(stderr, "done.\n");


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

    processMappedCs(sr, cs, mapc, mg, tabBases, tab_cpt, &nbErrors, &nbSearch, tabObs, mos, mobs, minCov, maxFracAlter);

    nbSamReadDone++;
    if(nbSamReadDone%100000 == 0){
      fprintf(stderr, "%d sam reads done.\n", nbSamReadDone);
    }
  }

  fprintf(stderr, "%ld errors in %ld bases.\n", nbErrors, nbSearch);

  fprintf(stderr, "Observations:\n");
  for(i=0 ; i<4 ; i++){
    fprintf(stderr, "%c\t%ld\n", tabBases[i], tabObs[tabBases[i]]);
  }

  FILE *nfObs;
  nfObs = fopen(argv[ARG_FILE_OBS], "w");
  mobs.print(nfObs);
  fclose(nfObs);

  return 0;
}



int processMappedCs(samRead &sr, CS_consensus &cs, circMapCovSimple &mapc, map<string, string> &mg, char *tabBases, char *tab_cpt, long *nbErrors, long *nbSearch, long *tabObs, mapOfStrand &mos, circMapObs &mobs, int minCoverage, double maxFracAlter){

  bool refined = false;

  string sCigar = sr.getCigar();
  if( sCigar.find("D")!=string::npos || sCigar.find("I")!=string::npos )
    return 0;

  //int minCoverage = 100;
  //double maxFracAlter = 0.01;
  int minScore = 100;
  //minScore = 70; // DEBOG !!!!

  int nbRef, nbTot;
  vector<int> vProbs;
  vector<int> vDepth;
  vector<int> vDivScores;

  string readID = sr.getQueryName();
  string chrID = sr.getRefName();

  //fprintf(stderr, "DEBOG: WORKING ON '%s'\n", readID.c_str());

  bool readMapsOnReverse = false;
  if( sr.readReverseStrand() == true ){
    cs.reverseCpt(tab_cpt);
    readMapsOnReverse = true;
    //return 0;
  }


  map<int,int> mc;
  // DEBOG
  //if( sr.getGenomicCoords(mc) == false )
  //return 0;
  if( sr.getMapOfCoords(mc)!= 0 )
    return 0;

  int nbMM = sr.getNbMM();
  if( nbMM > 1 ){ return -1; } // I do not even look at reads that have more than one mismatch
  if( nbMM==1 ){
    ;
    /*
    int offset = 0;
    int res = cs.refine(mc, mg[chrID], &offset);
    if( res != 0 ){
      return 0;
    }

    sr.setPos( sr.getPos() + offset );
    if( sr.getGenomicCoords(mc) == false )
      return 0;
    refined = true;
    */
  }


  cs.updateConsensus();
  string scs = cs.getConsensus();
  int csOK = cs.getCsProb(vProbs, vDepth, vDivScores);
  if(csOK<0){
    fprintf(stderr, "ABORTING THIS READ.\n");
    return 0; // DO SOMETHING!
  }

  int alnLg = sr.getSeq().size();

  map<int,int>::iterator itmc;
  for(itmc=mc.begin() ; itmc!=mc.end() ; itmc++){
    int posCs = itmc->first;

    // Exclude the first 5 and last 5 positions in the consensus
    // DEBOG: RE-DECOMMENT THIS PART AFTER FINISHED TEST FOR BASTIEN
    //if( posCs<5 || posCs>(alnLg-5) )
    //continue;
    
    int posGenomic = itmc->second;
    if( posGenomic<0) // = If I'm looking at a deletion
      continue;

    /////////////////////////////////////////////////////////////////////////////
    // Looking for polymorphic site
    if( mapc.isSuspect(chrID, posGenomic, minCoverage, maxFracAlter, tabBases, &nbRef, &nbTot)==true ){
      continue;
    }

    // Looking only at positions that have a high enough confidence level
    if( vProbs[posCs] < minScore )
      continue;

    // Set here the Minimum depth
    if( vDepth[posCs] < 3 )
      continue;

    // Looking only at places that have a perfect consensus (= no alternative base call whatsoever)
    if( vDivScores[posCs] > 0 ) // !!! DE-COMMENT THIS PART FOR PRODUCTION !!!
	continue;

    char cStrand = mos.m_[chrID].at(posGenomic);
    if( cStrand!='+' && cStrand!='-' )
      continue;

    bool posOnCpt = false;
    if( cStrand=='-' )
      posOnCpt = true;


    (*nbSearch)++;
    mobs.addObs(chrID, posGenomic);

    if(posGenomic<0){
      fprintf(stderr, "posGenomic: %d\nABORTING PROGRAM.\n", posGenomic);
      exit(1);
    }

    char refBase = mg[chrID].at(posGenomic);
    char csBase = scs[posCs];
    //char rnaBase = csBase;

    if( posOnCpt==true ){
      //rnaBase = tab_cpt[rnaBase];
      refBase = tab_cpt[refBase];
      csBase = tab_cpt[csBase];
      ;
    }
    tabObs[refBase]++;

    if( refBase != csBase ){ // && vDivScores[posCs]>78 ){
      // DEBOG: Looking specifically at cases with an ambiguous consensus
  
      fprintf(stdout, "nbRef: %d - nbTot: %d (%.4f)\n", nbRef, nbTot, (double)nbRef/(double)nbTot);
      //fprintf(stdout, "Refined: ");
      //if( refined==true ){ fprintf(stdout, "YES.\n"); } else { fprintf(stdout, "NO.\n"); }

      //char phasedCsBase = csBase;
      //if( posOnCpt==true )
      //phasedCsBase = tab_cpt[csBase];

      fprintf(stdout, "#%s\t%d\t%d\t%c\t%c\t%s\t%d\t%c\n", readID.c_str(), posCs, alnLg, refBase, csBase, chrID.c_str(), posGenomic, cStrand);
      string genomicSeq = sr.getGenomicSeq(mc, mg);
	
      // DEBOG: Add N nucleotides on both ends of the genomique sequence.
      int extraGenomic = 20;
      int minPosGenomic = mc.begin()->second;
      int maxPosGenomic = mc.rbegin()->second;

      int chrSize = mg[chrID].size();


      if( (minPosGenomic-extraGenomic) < 0){
	extraGenomic = minPosGenomic;
      }
      if( (maxPosGenomic+extraGenomic)>chrSize ){
	extraGenomic = (chrSize-maxPosGenomic) - 1;
      }



      string extendedGenomicSeq = mg[chrID].substr(minPosGenomic-extraGenomic, extraGenomic) + genomicSeq + mg[chrID].substr(maxPosGenomic+1, extraGenomic);
      string sTest = mg[chrID].substr(minPosGenomic-extraGenomic, alnLg+(2*extraGenomic));
	
      fprintf(stdout, "\n!!!!!!!!!!!!!\n");
      fprintf(stdout, 
	      "%s\n%s\n%*c%s\n",
	      sTest.c_str(),
	      extendedGenomicSeq.c_str(),
	      extraGenomic,
	      ' ',
	      scs.c_str()
	      );
      fprintf(stdout, "\n!!!!!!!!!!!!!!!\n");
	
      //fprintf(stderr, "genomicSeq:'%s'\n", genomicSeq.c_str());
      printDetailedMM(stdout, cs, genomicSeq, posCs);
      (*nbErrors)++;
      fflush(stdout);

    }// END of IF I'm looking at a valid mismatch

  }// END of going through all positions in mc

}



	/*
	// ADD HERE THE PART THAT TRIES TO REFINE THE ALIGNEMENT AND PRINTS THE REFINED ALIGNEMENT
	//if(queryName=="HSQ-7001360:223:H3L5LBCXX:1:2216:5049:8074")
	//fprintf(stderr, "break.\n");
	map<int, int> refinedMap;
	string refinedCs = "";
	string refinedGenomic = "";
	cs.refine(mc, posCs, mg[chrID], refinedMap, refinedCs);
	map<int,int>::iterator itref;
	for(itref=refinedMap.begin() ; itref!=refinedMap.end() ; itref++){
	  int pg = itref->second;
	  if(pg < 0){
	    refinedGenomic.push_back('-');
	  } else {
	    refinedGenomic.push_back(mg[chrID].at(pg));
	  }
	}
	fprintf(stdout, "+++++++++++++++++++++++++++\n");
	for(int ii=0 ; ii<refinedCs.size() ; ii++){
	  if( refinedCs[ii]==refinedGenomic[ii] ){
	    fprintf(stdout, "%c", ' ');
	  } else {
	    fprintf(stdout, "%c", '*');
	  }
	}
	fprintf(stdout, "\n");
	fprintf(stdout, "%s\n%s\n", refinedCs.c_str(), refinedGenomic.c_str());
	fprintf(stdout, "+++++++++++++++++++++++++++\n");
	fflush(stdout);
	*/




void printDetailedMM(FILE *nf, CS_consensus &cs, string &genomeSeq, int posCs){

  string scs = cs.getConsensus();

  for(int ii=0 ; ii<posCs ; ii++){
    fprintf(stdout, "%c", ' ');
  }
  fprintf(stdout, "%c\n", '*');
  
  fprintf(stdout, "%s\n%s\n", genomeSeq.c_str(), scs.c_str());

  fprintf(stdout, "\n-----\n");
  cs.prettyPrint(stdout, false);
  fprintf(stdout, "\n-----\n");
}
