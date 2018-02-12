/*

  This program that parses the sam output from the mapping of consensus sequences and calls candidate transcription errors (including indels).

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


using namespace std;

int updateIndelCovSimple(long *nbSearchIndel, vector<int> &vDepth, int minPos, int maxPos, int minRepeat, int maxRepeat);

int processMappedCs(samRead &sr, 
		    CS_consensus &cs, 
		    circMapCovSimple &mapc, 
		    C_circMapIndels &mapci,
		    long *nbIns, long *nbDel, long *nbSearchIndels,
		    map<string, string> &mg, 
		    char *tabBases, 
		    char *tab_cpt, 
		    long *nbErrors, long *nbSearch, long *tabObs, 
		    mapOfStrand &mos, 
		    circMapObs &mobs, 
		    int minCoverage, 
		    double maxFracAlter, 
		    int minRepeat, int maxRepeat, 
		    int nbToExcludeFromLeft, int nbToExcludeFromRight,
		    int minScore,
		    bool onlyPerfectCs,
		    map<string, GFFtranscript> &mt
		    );

int processIndel(
		  long *nbToIncrement, char indelType, 
		  //map<int, pair<int,int> >::iterator iti, 
		  int posInRead,
		  long gPos,
		  int minPos, int maxPos, 
		  vector<int> &vDepth, int minRepeat, int maxRepeat, 
		  string &chrID, 
		  C_circMapIndels &mapci, 
		  int minCov, double maxFracAlter,
		  map<string, GFFtranscript> &mt,
		  char *tab_cpt
		 );

int processIndels(
		  std::map<string, string> &mg,
		  mapOfStrand &mos,
		  samRead &sr,
		  C_circMapIndels &mapci,
		  long *nbIns, long *nbDel, long *nbSearchIndels,
		  int minCoverage, double maxFracAlter,
		  int minRepeat, int maxRepeat,
		  int nbToExcludeFromLeft, int nbToExcludeFromRight,
		  vector<int> &vDepth,
		  map<string, GFFtranscript> &mt,
		  char *tab_cpt
		  );



void printDetailedMM(FILE *nf, CS_consensus &cs, string &genomeSeq, int posCs);


//#define ARG_FILE_DEBREAKED 1
#define ARG_FILE_BLAT_POS 1
#define ARG_FILE_CS_POS 2
#define ARG_FILE_SAM 3
#define ARG_KEEP_ONLY_UNIQUE 4

#define ARG_FILE_GENOME 5
#define ARG_FILE_GFF 6
#define ARG_FILE_GFF_TYPE 7
#define ARG_FILE_GFF_KEY 8

#define ARG_FILE_MAP 9
#define ARG_FILE_MAP_INDELS 10

#define ARG_MIN_COV 11
#define ARG_MAX_FRAC_ALTER 12
#define ARG_FILE_OBS 13

#define ARG_MIN_REPEAT 14
#define ARG_MAX_REPEAT 15

#define ARG_NB_EXCLUDE_LEFT 16
#define ARG_NB_EXCLUDE_RIGHT 17

#define ARG_MIN_SCORE 18

#define ARG_ONLY_PERFECT_CS 19

#define ARG_CHRID 20

#define ARG_FILE_READS 21
#define ARG_FILE_MATES 22



int main(int argc, char *argv[]){

  if(argc<2){
    fprintf(stderr, "%s: calls transcription errors from circseq data.\n", argv[0]);
    fprintf(stderr, "Usage:\n%s blat-mapping.pos lib-CS.pos map.sam keepOnlyUniqueMapped [T/F] ref-genome.fa ref.gtf [GFF/GTF] gffKey map.tab map-indels.tab minCov maxFracAlter file-obs minRepeat maxRepeat nbExcludeLEft nbExlcudeRight minScore onlyPerfectCs [T/F] chrID reads.fastq [mates.fastq]\n", argv[0]);
    return 1;
  }

  int i;

  long nbErrors = 0;
  long nbSearch = 0;

  long nbIns = 0;
  long nbDel = 0;
  long nbSearchIndels = 0;

  int minCov = atoi(argv[ARG_MIN_COV]);
  double maxFracAlter = atof(argv[ARG_MAX_FRAC_ALTER]);

  int minRepeat = atoi(argv[ARG_MIN_REPEAT]);
  int maxRepeat = atoi(argv[ARG_MAX_REPEAT]);

  int nbToExcludeFromLeft = atoi(argv[ARG_NB_EXCLUDE_LEFT]);
  int nbToExcludeFromRight = atoi(argv[ARG_NB_EXCLUDE_RIGHT]);
  int minScore = atoi(argv[ARG_MIN_SCORE]);

  char ctmp = argv[ARG_ONLY_PERFECT_CS][0];
  bool onlyPerfectCs = false;
  if(ctmp=='T'){
    onlyPerfectCs = true;
  }

  string chrOnly(argv[ARG_CHRID]);

  char tabBases[4] = {'A', 'T', 'C', 'G'};
  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);

  long tabObs['z'+'Z'];
  for(i=0 ; i<'z'+'Z' ; i++){
    tabObs[i] = 0;
  }

  bool keepOnlyUniqueMapped = false;
  if( argv[ARG_KEEP_ONLY_UNIQUE][0]=='T' ){
    keepOnlyUniqueMapped = true;
  }

  map<string, string> mg;
  fprintf(stderr, "Loading genome in memory (from '%s') ...", argv[ARG_FILE_GENOME]);
  loadFasta(mg, argv[ARG_FILE_GENOME], chrOnly);
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

  vector<string> vSeqNames;
  vSeqNames.push_back(chrOnly);  

  getAllTranscripts(argv[ARG_FILE_GFF], 0, mt, gffType, gffKey, 1, mg, tab_cpt, vSources, vFeatures, vSeqNames);
  mapOfStrand mos(mg, mt);


  fprintf(stderr, "Initializating map of coverage ...");
  circMapCovSimple mapc(mg);
  fprintf(stderr, " done.\n");
  fprintf(stderr, "Loading values into map of coverage (from '%s') ...", argv[ARG_FILE_MAP]);
  mapc.load(argv[ARG_FILE_MAP], chrOnly);
  fprintf(stderr, "done.\n");

  fprintf(stderr, "Initializating the map of indels ...");
  fflush(stderr);
  C_circMapIndels mapci(mg);
  fprintf(stderr, " done.\n");
  fprintf(stderr, "Loading values into map of indels (from '%s') ...", argv[ARG_FILE_MAP_INDELS]);
  mapci.load(argv[ARG_FILE_MAP_INDELS], chrOnly);
  fprintf(stderr, " done.\n");


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

    //Only reads mapped to the chromosome being studied (for memory optimization)
    if( sr.getRefName() != chrOnly )
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

    processMappedCs(
		    sr, 
		    cs, 
		    mapc, 
		    mapci,
		    &nbIns, &nbDel, &nbSearchIndels,
		    mg, 
		    tabBases, 
		    tab_cpt, 
		    &nbErrors, &nbSearch, 
		    tabObs, mos, mobs, 
		    minCov, maxFracAlter, 
		    minRepeat, maxRepeat, 
		    nbToExcludeFromLeft, nbToExcludeFromRight, 
		    minScore,
		    onlyPerfectCs,
		    mt
		    );

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

  fprintf(stderr, "%ld insertion / %ld deletions / %ld indelObservations\n", nbIns, nbDel, nbSearchIndels);

  FILE *nfObs;
  nfObs = fopen(argv[ARG_FILE_OBS], "w");
  mobs.print(nfObs);
  fclose(nfObs);

  return 0;
}



int processMappedCs(samRead &sr, 
		    CS_consensus &cs, 
		    circMapCovSimple &mapc, 
		    C_circMapIndels &mapci,
		    long *nbIns, long *nbDel, long *nbSearchIndels,
		    map<string, string> &mg, 
		    char *tabBases, 
		    char *tab_cpt, 
		    long *nbErrors, long *nbSearch, long *tabObs, 
		    mapOfStrand &mos, 
		    circMapObs &mobs, 
		    int minCoverage, 
		    double maxFracAlter, 
		    int minRepeat, int maxRepeat, 
		    int nbToExcludeFromLeft, int nbToExcludeFromRight,
		    int minScore,
		    bool onlyPerfectCs,
		    map<string, GFFtranscript> &mt		   
		    ){

  //bool refined = false;

  string sCigar = sr.getCigar();

  //int minCoverage = 100;
  //double maxFracAlter = 0.01;
  //int minScore = 100;
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

  double csId = cs.getPrctId();
  if( csId<0.985 ){
    return 0;
  }

  // DEBOG !
  //cs.prettyPrint(stderr, false);
  //fprintf(stderr, "Consensus identity prct: %.4f\n", cs.getPrctId());

  map<int,int> mc;
  if( sr.getMapOfCoords(mc) != 0 )
    return 0;

  int eDist = sr.getEditDistance();
  if( eDist > 1 ){ return -1; } // I do not even look at reads that have more than one mismatch or one mismatch + one indel


  cs.updateConsensus();
  string scs = cs.getConsensus();
  int csOK = cs.getCsProb(vProbs, vDepth, vDivScores);
  if(csOK<0){
    fprintf(stderr, "ABORTING THIS READ.\n");
    return 0; // DO SOMETHING!
  }

  int alnLg = sr.getSeq().size();


  // I can keep this part to prevent calling mismatches from indel-containing reads
  if( sCigar.find("D")!=string::npos || sCigar.find("I")!=string::npos ){
    // Here, process indels.
    long iNull;
    processIndels(
		  mg,
		  mos,
		  sr, 
		  mapci, 
		  nbIns, nbDel, &iNull,  // replace &iNull with nbSearchIndels if I do not pass indel-containing reads thru the normal part of the pipeline.
		  minCoverage, maxFracAlter, 
		  minRepeat, maxRepeat, 
		  nbToExcludeFromLeft, nbToExcludeFromRight,
		  vDepth,
		  mt,
		  tab_cpt
		  );
    //return 0; // <- de-comment this part to prevent searching for base-subs in reads that contain indels.
  }


  map<int,int>::iterator itmc;
  for(itmc=mc.begin() ; itmc!=mc.end() ; itmc++){
    int posCs = itmc->first;

    // Exclude the first 5 and last 5 positions in the consensus
    if( posCs<nbToExcludeFromLeft || posCs>((alnLg-nbToExcludeFromRight)-1) )
      continue;
    
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
    if( vDepth[posCs]<minRepeat || vDepth[posCs]>maxRepeat )
      continue;

    char cStrand = mos.m_[chrID].at(posGenomic);
    //if( cStrand!='+' && cStrand!='-' )
    //continue;

    // Looking only at places that have a perfect consensus (= no alternative base call whatsoever) if this option is set to true
    if( onlyPerfectCs==true && vDivScores[posCs] > 0 )
      continue;

    if( mapci.isSuspect(chrID, posGenomic, minCoverage, maxFracAlter)==false ){
      (*nbSearchIndels)++;
    }


    bool posOnCpt = false;
    if( cStrand=='-' )
      posOnCpt = true;


    (*nbSearch)++;
    if( posGenomic<0 || posGenomic>mg[chrID].size() ){
      fprintf(stderr, " ERROR - POS GENOMIC OUT OF RANGE.\n");
      return 0;
    }

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
  
      //fprintf(stdout, "nbRef: %d - nbTot: %d (%.4f)\n", nbRef, nbTot, (double)nbRef/(double)nbTot);
      //fprintf(stdout, "Refined: ");
      //if( refined==true ){ fprintf(stdout, "YES.\n"); } else { fprintf(stdout, "NO.\n"); }

      //char phasedCsBase = csBase;
      //if( posOnCpt==true )
      //phasedCsBase = tab_cpt[csBase];

      string gID = "NO_GENE";
      vector<GFFtranscript> vt;
      getAllTranscriptsAtPosition(vt, mt, chrID, posGenomic);
      if( vt.size() > 0 ){
	gID = vt[0].id_;
      }

      fprintf(stdout, 
	      "#BASE_SUB\t%s\t%d\t%d\t%c\t%c\t%s\t%d\t%c\t%s\n", 
	      readID.c_str(), posCs, alnLg, refBase, csBase, chrID.c_str(), posGenomic, cStrand, gID.c_str());

      /* ********************************************************************************************************************************* */
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
      /* ********************************************************************************************************************************* */

    }// END of IF I'm looking at a valid mismatch

  }// END of going through all positions in mc

}



int processIndels(
		  map<string, string> &mg,
		  mapOfStrand &mos,
		  samRead &sr,
		  C_circMapIndels &mapci,
		  long *nbIns,
		  long *nbDel,
		  long *nbSearchIndels,
		  int minCoverage, double maxFracAlter,
		  int minRepeat, int maxRepeat,
		  int nbToExcludeFromLeft, int nbToExcludeFromRight,
		  vector<int> &vDepth,
		  map<string, GFFtranscript> &mt,
		  char *tab_cpt
		  ){

  int cov, nbInsAtPos, nbDelAtPos;
  char bFrom, bTo, strand;
  string sFrom, sTo;
  int indelSize;
  double rr;
  map<int, pair<int,int> > mIns;
  map<int, pair<int,int> > mDel;

  int posInRead, gPos, size;
  pair<int,int> pp;

  sr.getIndels(mIns, mDel);

  string chrID = sr.getRefName();
  int readLg = sr.getSeqLg();
  string readSeq = sr.getSeq();

  int minPos = nbToExcludeFromLeft;
  int maxPos = (readLg - nbToExcludeFromRight) - 1;

  map<int, pair<int,int> >::iterator iti;
  for(iti=mIns.begin() ; iti!=mIns.end() ; iti++){
    posInRead = iti->first;
    gPos = (iti->second).first;
    indelSize = (iti->second).second;
    int foundIns = processIndel(nbIns, 'I', posInRead, gPos, minPos, maxPos, vDepth, minRepeat, maxRepeat, chrID, mapci, minCoverage, maxFracAlter, mt, tab_cpt);
    strand = mos.m_[chrID].at(gPos);
    if( foundIns>0 ){
      bFrom = '-';
      bTo = readSeq[posInRead];
      sFrom = string(indelSize, '-');
      sTo = readSeq.substr(posInRead, indelSize);

      if(strand=='-'){
	sTo = my_reverse_cpt(sTo, tab_cpt);
      }

      string gID = "NO_GENE";
      vector<GFFtranscript> vt;
      getAllTranscriptsAtPosition(vt, mt, chrID, gPos);
      if( vt.size() > 0 ){
	gID = vt[0].id_;
      }

      fprintf(stdout, 
	      "#INSERTION\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%c\t%s\n", 
	      sr.getQueryName().c_str(), 
	      posInRead,
	      readLg,
	      sFrom.c_str(),
	      sTo.c_str(),
	      chrID.c_str(),
	      gPos,
	      strand,
	      gID.c_str()
	      );
      updateIndelCovSimple(nbSearchIndels, vDepth, minPos, maxPos, minRepeat, maxRepeat);
    }
  }

  map<int, pair<int,int> >::iterator itd;
  for(itd=mDel.begin() ; itd!=mDel.end() ; itd++){
    posInRead = itd->first;
    gPos = (itd->second).first;
    indelSize = (itd->second).second;
    strand = mos.m_[chrID].at(gPos);
    int foundDel = processIndel(nbDel, 'D', posInRead, gPos, minPos, maxPos, vDepth, minRepeat, maxRepeat, chrID, mapci, minCoverage, maxFracAlter, mt, tab_cpt);
    if( foundDel>0 ){
      sFrom = mg[chrID].substr(gPos+1, indelSize);
      sTo = string(indelSize, '-');

      if(strand=='-'){
	sFrom = my_reverse_cpt(sFrom, tab_cpt);
      }

      string gID = "NO_GENE";
      vector<GFFtranscript> vt;
      getAllTranscriptsAtPosition(vt, mt, chrID, gPos);
      if( vt.size() > 0 ){
	gID = vt[0].id_;
      }

      //fprintf(stdout, "#DELETION\t%s\n", sr.getQueryName().c_str());
      fprintf(stdout,
	      "#DELETION\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%c\t%s\n", 
	      sr.getQueryName().c_str(), 
	      posInRead,
	      readLg,
	      sFrom.c_str(),
	      sTo.c_str(),
	      chrID.c_str(),
	      gPos,
	      strand,
	      gID.c_str()
	      );

      updateIndelCovSimple(nbSearchIndels, vDepth, minPos, maxPos, minRepeat, maxRepeat);
    }
  }



  return 0;
}


int updateIndelCovSimple(long *nbSearchIndels, vector<int> &vDepth, int minPos, int maxPos, int minRepeat, int maxRepeat){

  int i;
  int nb = vDepth.size();

  for(i=0 ; i<nb ; i++){
    int dd = vDepth[i];
    if(i>=minPos && i<=maxPos && dd>=minRepeat && dd<=maxRepeat){
      (*nbSearchIndels)++;
    }
  }

}


int processIndel(
		  long *nbToIncrement, char indelType, 
		  //map<int, pair<int,int> >::iterator iti, 
		  int posInRead,
		  long gPos,
		  int minPos, int maxPos, 
		  vector<int> &vDepth, int minRepeat, int maxRepeat, 
		  string &chrID, 
		  C_circMapIndels &mapci, 
		  int minCov, double maxFracAlter,
		  map<string, GFFtranscript> &mt,
		  char *tab_cpt
		 ){

  long cov, nbInsAtPos, nbDelAtPos;
  double rr;

  //int posInRead = iti->first;
  if( posInRead<minPos || posInRead>maxPos )
    return 0;

  if( vDepth[posInRead]<minRepeat || vDepth[posInRead]>maxRepeat )
    return 0;

  // Special indels: check the coverage of the prvious and next positions?
  int readLg = vDepth.size();
  if( posInRead<1 || posInRead>(readLg-1) || vDepth[posInRead-1]<minRepeat || vDepth[posInRead+1]<minRepeat )
    return 0;

  //pair<int, int> pp = iti->second;
  //gPos = pp.first;
  if( mapci.isSuspect(chrID, gPos, minCov, maxFracAlter) ){
    return 0;
  }

  (*nbToIncrement)++;

  return 1;
}


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
