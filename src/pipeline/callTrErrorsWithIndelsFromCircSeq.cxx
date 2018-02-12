/*

  This program that parses the sam output from the mapping of consensus sequences and calls candidate transcription errors (including indels).

 */

//#define _DEBOG_ 1

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

int updateIndelCovSimple(long *nbSearchIndel, vector<int> &vDepth, int minPos, int maxPos, int minRepeat, int maxRepeat);

int processMappedCs(samRead &sr, 
		    CS_consensus &cs, 
		    circMapGenome & genomeMap,
		    long *nbIns, long *nbDel, long *nbSearchIndels,
		    map<string, string> &mg, 
		    char *tabBases, 
		    char *tab_cpt, 
		    long *nbErrors, long *nbSearch, long *tabObs, 
		    mapOfStrand &mos, 
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
		  circMapGenome & genomeMap,
		  map<string, map<long,circMapGenomePos> >::iterator itc,
		  int minCov, double maxFracAlter,
		  map<string, GFFtranscript> &mt,
		  char *tab_cpt
		 );

int processIndels(
		  std::map<string, string> &mg,
		  mapOfStrand &mos,
		  samRead &sr,
		  circMapGenome & genomeMap,
		  map<string, map<long,circMapGenomePos> >::iterator itc,
		  long *nbIns, long *nbDel, long *nbSearchIndels,
		  int minCoverage, double maxFracAlter,
		  int minRepeat, int maxRepeat,
		  int nbToExcludeFromLeft, int nbToExcludeFromRight,
		  vector<int> &vDepth,
		  map<string, GFFtranscript> &mt,
		  char *tab_cpt
		  );



bool positionIsCallable(
			map<int,int>::iterator &itmc, 
			int nbToExcludeFromLeft, int nbToExcludeFromRight, int alnLg,
			circMapGenome &genomeMap,
			map<string, map<long,circMapGenomePos> >::iterator itc,
			int minCoverage,
			double maxFracAlter,
			vector<int> &vProbs,
			int minScore,
			vector<int> &vDivScores,
			int minRepeat, int maxRepeat,
			vector<int> &vDepth,
			bool onlyPerfectCs,
			char *tabBases
			);


void updateOK(string &scs, vector<bool> &vok);

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

#define ARG_FILE_READS 20
#define ARG_FILE_MATES 21



int main(int argc, char *argv[]){

  if(argc<2){
    fprintf(stderr, "%s: calls transcription errors from circseq data.\n", argv[0]);
    fprintf(stderr, "Usage:\n%s blat-mapping.pos lib-CS.pos map.sam keepOnlyUniqueMapped [T/F] ref-genome.fa ref.gtf [GFF/GTF] gffKey map.tab map-indels.tab minCov maxFracAlter file-obs minRepeat maxRepeat nbExcludeLEft nbExlcudeRight minScore onlyPerfectCs [T/F] reads.fastq [mates.fastq]\n", argv[0]);
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
  loadFasta(mg, argv[ARG_FILE_GENOME]);
  fprintf(stderr, " done.\n");

#ifdef _DEBOG_
  string seqTmp = "";
  for(int iip=4435138-10 ; iip<4435138+10 ; iip++){
    char baseTmp = mg["V"][iip];
    fprintf(stdout, "%d:%c\n", iip, baseTmp);
  }
  fprintf(stdout, "%s\n", seqTmp.c_str());
#endif

  circMapGenome genomeMap;

  string gffType(argv[ARG_FILE_GFF_TYPE]);
  string gffKey(argv[ARG_FILE_GFF_KEY]);

  map<string, GFFtranscript> mt;
  vector<string> vSources;
  vector<string> vFeatures;
  vector<string> vChr;
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

  getAllTranscripts(argv[ARG_FILE_GFF], 0, mt, gffType, gffKey, 1, mg, tab_cpt, vSources, vFeatures, vChr);
  mapOfStrand mos(mg, mt);


  fprintf(stderr, "Loading values into map of coverage (from '%s') ...", argv[ARG_FILE_MAP]);
  genomeMap.loadMapCoverage(argv[ARG_FILE_MAP], false);
  fprintf(stderr, "done.\n");

  fprintf(stderr, "Loading values into map of indels (from '%s') ...", argv[ARG_FILE_MAP_INDELS]);
  genomeMap.loadMapIndels(argv[ARG_FILE_MAP_INDELS]);
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

    // Keep only uniquely matched reads
    int iFlag = sr.getFlag();
    if( iFlag != 0 && iFlag != 16 )
      continue;

    if( keepOnlyUniqueMapped==true && sr.mOpt_["NH"]!="1" )
      continue;

    // Skip mapped reads that have more than 3 components to their CIGAR string (for example 40M1I1N20M would be skipped).
    // This should allow skipping all these strange cases where tophat calls an indels while there is a perfect match...
    int cigarDepth = sr.getCigarDepth();
    if( cigarDepth > 3 )
      continue;

    ///////////////////////////////////////////////////////////////////
    //                DEBOG ONLY - REMOVE FOR PRODUCTION 
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //if( sr.isOnIntron()==true )
    //continue;

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

    cs.setLimitHomo();
    //cs.setLimitHomo(0, 0);

    processMappedCs(
		    sr, 
		    cs, 
		    genomeMap,
		    &nbIns, &nbDel, &nbSearchIndels,
		    mg, 
		    tabBases, 
		    tab_cpt, 
		    &nbErrors, &nbSearch, 
		    tabObs, mos, 
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
  genomeMap.printObservations(nfObs);
  fclose(nfObs);

  return 0;
}



int processMappedCs(samRead &sr, 
		    CS_consensus &cs, 
		    circMapGenome & genomeMap,
		    long *nbIns, long *nbDel, long *nbSearchIndels,
		    map<string, string> &mg, 
		    char *tabBases, 
		    char *tab_cpt, 
		    long *nbErrors, long *nbSearch, long *tabObs, 
		    mapOfStrand &mos, 
		    int minCoverage, 
		    double maxFracAlter, 
		    int minRepeat, int maxRepeat, 
		    int nbToExcludeFromLeft, int nbToExcludeFromRight,
		    int minScore,
		    bool onlyPerfectCs,
		    map<string, GFFtranscript> &mt		   
		    ){

  //bool refined = false;

  int homoLeft, homoRight;
  cs.getHomoSizes(&homoLeft, &homoRight);

  if(homoLeft > nbToExcludeFromLeft){ nbToExcludeFromLeft = homoLeft; }
  if(homoRight > nbToExcludeFromLeft){ nbToExcludeFromRight = homoRight; }


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

  map<string, map<long,circMapGenomePos> >::iterator itc = genomeMap.mg_.find(chrID);
  if( itc == genomeMap.mg_.end() ){
    fprintf(stderr, "WARNING: '%s' not found in genome map. This should never happen!\n", chrID.c_str());
    return -1;
  }

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
		  genomeMap, 
		  itc,
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


  vector<bool> vok(scs.size(), true);
  map<int,int>::iterator itmc;
  int posCs, posGenomic;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(itmc=mc.begin() ; itmc!=mc.end() ; itmc++){
    posCs = itmc->first;
    vok[posCs] = positionIsCallable(itmc, 
				    nbToExcludeFromLeft, nbToExcludeFromRight, alnLg,
				    genomeMap,
				    itc,
				    minCoverage,
				    maxFracAlter,
				    vProbs,
				    minScore,
				    vDivScores,
				    minRepeat, maxRepeat, vDepth,
				    onlyPerfectCs,
				    tabBases
				    );
  }
  updateOK(scs, vok);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////


  for(itmc=mc.begin() ; itmc!=mc.end() ; itmc++){
    posCs = itmc->first;
    posGenomic = itmc->second;


    // Exclude the first 5 and last 5 positions in the consensus
    //if( posCs<nbToExcludeFromLeft || posCs>((alnLg-nbToExcludeFromRight)-1) )
    //continue;
    
    //int posGenomic = itmc->second;
    //if( posGenomic<0) // = If I'm looking at a deletion
	//continue;


    /////////////////////////////////////////////////////////////////////////////
    // Looking for polymorphic site
    //if( genomeMap.isSuspect_cov(itc, posGenomic, minCoverage, maxFracAlter, tabBases, &nbRef, &nbTot)==true ){
    //continue;
    //}

    // Looking only at positions that have a high enough confidence level
    //if( vProbs[posCs] < minScore )
    //continue;

    // Set here the Minimum depth
    //if( vDepth[posCs]<minRepeat || vDepth[posCs]>maxRepeat )
    //continue;

    //if( cStrand!='+' && cStrand!='-' )
    //continue;

    // Looking only at places that have a perfect consensus (= no alternative base call whatsoever) if this option is set to true
    //if( onlyPerfectCs==true && vDivScores[posCs] > 0 )
    //continue;

    //if( genomeMap.isSuspect_indels(itc, posGenomic, minCoverage, maxFracAlter)==false ){
    //(*nbSearchIndels)++;
    //}

    if( vok[posCs] == false )
      continue;

    (*nbSearchIndels)++;

    char cStrand = mos.m_[chrID].at(posGenomic);

    bool posOnCpt = false;
    if( cStrand=='-' )
      posOnCpt = true;

#ifdef _DEBOG_
    fprintf(stdout, "[%d -> %d] [%c] [%c:%c]\n", posCs, posGenomic, cStrand, scs[posCs], mg[chrID].at(posGenomic));
#endif

    (*nbSearch)++;
    if( posGenomic<0 || posGenomic>mg[chrID].size() ){
      fprintf(stderr, " ERROR - POS GENOMIC OUT OF RANGE.\n");
      return 0;
    }

    genomeMap.addObservation(itc, posGenomic);

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
		  circMapGenome &genomeMap,
		  map<string, map<long,circMapGenomePos> >::iterator itc,
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
    int foundIns = processIndel(nbIns, 'I', posInRead, gPos, minPos, maxPos, vDepth, minRepeat, maxRepeat, chrID, genomeMap, itc, minCoverage, maxFracAlter, mt, tab_cpt);
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
    int foundDel = processIndel(nbDel, 'D', posInRead, gPos, minPos, maxPos, vDepth, minRepeat, maxRepeat, chrID, genomeMap, itc, minCoverage, maxFracAlter, mt, tab_cpt);
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
		  long *nbToIncrement, 
		  char indelType, 
		  //map<int, pair<int,int> >::iterator iti, 
		  int posInRead,
		  long gPos,
		  int minPos, int maxPos, 
		  vector<int> &vDepth, int minRepeat, int maxRepeat, 
		  string &chrID, 
		  circMapGenome &genomeMap,
		  map<string, map<long,circMapGenomePos> >::iterator itc,
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
  if( genomeMap.isSuspect_indels(itc, gPos, minCov, maxFracAlter) ){
    return 0;
  }
  // For insertions, I also look at the position -1
  if( indelType=='I' && genomeMap.isSuspect_indels(itc, gPos-1, minCov, maxFracAlter) ){
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


bool positionIsCallable(
			map<int,int>::iterator &itmc, 
			int nbToExcludeFromLeft, int nbToExcludeFromRight, int alnLg,
			circMapGenome &genomeMap,
			map<string, map<long,circMapGenomePos> >::iterator itc,
			int minCoverage,
			double maxFracAlter,
			vector<int> &vProbs,
			int minScore,
			vector<int> &vDivScores,
			int minRepeat, int maxRepeat,
			vector<int> &vDepth,
			bool onlyPerfectCs,
			char *tabBases
			){

  int posCs = itmc->first;
  int posGenomic = itmc->second;
  int nbRef, nbTot;

  if( posGenomic<0) // = If I'm looking at a deletion
    return false;    

  if( posCs<nbToExcludeFromLeft || posCs>((alnLg-nbToExcludeFromRight)-1) )
    return false;

  /////////////////////////////////////////////////////////////////////////////
  // Looking for polymorphic site
  if( genomeMap.isSuspect_cov(itc, posGenomic, minCoverage, maxFracAlter, tabBases, &nbRef, &nbTot)==true ){
    return false;
  }

  // Looking only at positions that have a high enough confidence level
  if( vProbs[posCs] < minScore )
    return false;

  // Set here the Minimum depth
  if( vDepth[posCs]<minRepeat || vDepth[posCs]>maxRepeat )
    return false;


  // Looking only at places that have a perfect consensus (= no alternative base call whatsoever) if this option is set to true
  if( onlyPerfectCs==true && vDivScores[posCs] > 0 )
    return false;

  return true;
}


void updateOK(string &scs, vector<bool> &vok){
  int csPos;
  int lg = scs.size();

  vok[0] = false;
  for(csPos=1 ; csPos<lg && scs[csPos]==scs[0] ; csPos++){
    vok[csPos] = false;
  }

  if(csPos>=lg)
    return;

  int startCurrentH = csPos;
  bool currentOK = vok[csPos];
  for(csPos=startCurrentH+1 ; csPos<lg ; csPos++){
    if(scs[csPos]==scs[startCurrentH]){
      if(vok[csPos]==false){ currentOK=false; }
    } else {
      for(int i=startCurrentH ; i<csPos ; i++){
	vok[i] = currentOK;
      }

      if(csPos<lg){
	startCurrentH = csPos;
	currentOK = vok[csPos];
      }

    }
  }

  vok[lg-1] = false;
  for(csPos=lg-2 ; csPos>=0 && scs[csPos]==scs[lg-1] ; csPos--){
    vok[csPos] = false;
  }

}
