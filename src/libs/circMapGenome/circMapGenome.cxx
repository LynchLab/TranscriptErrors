#ifndef _CIRCMAPGENOME_CXX_
#define _CIRCMAPGENOME_CXX_

#include "circMapGenome.hxx"

#include <stdlib.h>

void printMap(FILE *nf, map<int,int> &m);

//////////////////////////////////
//   circMapGenomePos functions    //
//////////////////////////////////


circMapGenomePos::circMapGenomePos(void){
  this->init();
}


circMapGenomePos::circMapGenomePos(simplePosCov &spc){
  this->init(spc);
}

circMapGenomePos::circMapGenomePos(C_indelPos &ip){
  this->init(ip);
}

circMapGenomePos::circMapGenomePos(long nbObs){
  this->init(nbObs);
}


char circMapGenomePos::initRefBase(map<string,string> &mg, string chrID, long gPos, mapOfStrand &mos, char *tab_cpt){

  char refBase = mg[chrID].at(gPos);
  char cStrand = mos.m_[chrID].at(gPos);
  if( cStrand!='+' && cStrand!='-' ){
    refBase = '?';
  } else {
    if( cStrand=='-' ){
      refBase = tab_cpt[refBase];
    }
  }
  this->spc_.setRefBase(refBase);
  return refBase;
}

void circMapGenomePos::init(void){
  simplePosCov spc;
  C_indelPos ip;
  this->spc_ = spc;
  this->ip_ = ip;
  this->nbObs_ = 0;
}


void circMapGenomePos::init(simplePosCov &spc){
  C_indelPos ip;
  this->spc_ = spc;
  this->ip_ = ip;
  this->nbObs_ = 0;
}

void circMapGenomePos::init(C_indelPos &ip){
  simplePosCov spc;
  this->spc_ = spc;
  this->ip_ = ip;
  this->nbObs_ = 0;
}

void circMapGenomePos::circMapGenomePos::init(long nbObs){
  simplePosCov spc;
  C_indelPos ip;
  this->spc_ = spc;
  this->ip_ = ip;
  this->nbObs_ = nbObs;
}




void circMapGenomePos::set(simplePosCov &spc){
  this->spc_ = spc;
}

void circMapGenomePos::set(C_indelPos &ip){
  this->ip_ = ip;
}

void circMapGenomePos::set(long nbObs){
  this->nbObs_ = nbObs;
}

void circMapGenomePos::add(C_indelPos &ip){
  this->ip_.add(ip);
}

void circMapGenomePos::add(simplePosCov &spc){
  this->spc_.add(spc);
}

void circMapGenomePos::add(long nbObs){
  this->nbObs_ += nbObs;
}






void circMapGenomePos::setRefBase(char refBase){
  this->spc_.setRefBase(refBase);
}


//???////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//                                 circMapGenome functions                                 //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////


circMapGenome::circMapGenome(void){
  this->init();
}


void circMapGenome::init(void){
  map<string, map<long,circMapGenomePos> > m;
  this->mg_ = m;
}


// BEGIN FUNCTIONS THAT DEAL WITH COVERGE


void circMapGenome::loadMapCoverage(char *fname, bool add){
  string line, stmp, chrID;
  int i;
  int tabCor['z'+'Z'];

  ifstream fic(fname);
  getline(fic, line);
  istringstream issh(line);

  getline(issh, stmp, '\t'); // "chrID"
  getline(issh, stmp, '\t'); // "pos"
  getline(issh, stmp, '\t'); // "refBase"

  for(i=0 ; i<4 ; i++){
    getline(issh, stmp, '\t');
    tabCor[tolower(stmp[0])] = i;
    tabCor[toupper(stmp[0])] = i;
  }

  map<string, map<long,circMapGenomePos> >::iterator itc;
  map<long,circMapGenomePos>::iterator itp;

  string prevChrID = "";

  while( getline(fic, line) ){
    istringstream iss(line);
    getline(iss, chrID, '\t');
    getline(iss, stmp, '\t');
    long gPos = atol(stmp.c_str());

    getline(iss, stmp, '\t');
    char refBase = stmp[0];

    if(chrID!=prevChrID){
      itc = this->mg_.find(chrID);
      if(itc==this->mg_.end()){
	map<long,circMapGenomePos> mc;
	itc = (this->mg_.insert( map<string, map<long,circMapGenomePos> >::value_type(chrID, mc) )).first;
      }
      prevChrID = chrID;
    }

    simplePosCov spc(refBase);
    for(i=0 ; i<4 ; i++){
      getline(iss, stmp, '\t');
      int nb = atoi(stmp.c_str());
      spc.set(i, nb);
      //this->mc_[chrID].at(gPos).set(i, nb);
    }

    itp = (itc->second).find(gPos);
    if( itp==(itc->second).end() ){
      circMapGenomePos gp;
      gp.set(spc);
      (itc->second).insert( map<long,circMapGenomePos>::value_type(gPos, gp) );
    } else {
      if( add==true ){
	(itp->second).add(spc);
      } else {
	(itp->second).set(spc);
      }
    }

  }

}


void circMapGenome::printMapCoverage(FILE *nf, const char *tabBases){
  int i;
  fprintf(nf, "%s\t%s\t%s", "scafID", "pos", "refBase");
  for(i=0 ; i<4 ; i++){
    fprintf(nf, "\t%c", tabBases[i]);
  }
  fprintf(nf, "\n");

  map<string, map<long,circMapGenomePos> >::const_iterator itc;
  for(itc=this->mg_.begin() ; itc!=this->mg_.end() ; itc++){
    string chrID = itc->first;
    map<long,circMapGenomePos> mp = itc->second;

    map<long,circMapGenomePos>::iterator itp;
    for(itp=mp.begin() ; itp!=mp.end() ; itp++){
      long gPos = itp->first;
      simplePosCov sc = (itp->second).spc_;
      fprintf(nf, "%s\t%d\t", chrID.c_str(), gPos);
      sc.print(nf);
      fprintf(nf, "\n");
    }
  }
}


int circMapGenome::updateMapCoverage(const char *fileSam, const int *tabCor, const int minScore){
  samRead sr;
  ifstream ficSam(fileSam);
  while( sr.readNext(ficSam) ){
    this->processSamRead_forCoverage(sr, tabCor, minScore);
  }
}




int circMapGenome::processSamRead_forCoverage(const samRead &sr, const int* tabCor, const int minScore){
  string chrID = sr.getRefName();
  map<string, map<long,circMapGenomePos> >::iterator itc = this->mg_.find(chrID);
  if( itc==this->mg_.end() ){
    map<long,circMapGenomePos> mtmp;
    itc = (this->mg_.insert( map<string, map<long,circMapGenomePos> >::value_type(chrID, mtmp) )).first;
  }
  this->processSamRead_forCoverage(sr, itc, tabCor, minScore);
}


int circMapGenome::processSamRead_forCoverage(
					      const samRead &sr, 
					      map<string,map<long,circMapGenomePos> >::iterator itc, 
					      const int* tabCor, 
					      const int minScore
					      ){
  
  string readSeq = sr.getSeq();
  string readID = sr.getQueryName();
  string chrID = sr.getRefName();

  map<int,int> mc;
  //sr.getGenomicCoords(mc);
  sr.getMapOfCoords(mc);

  //fprintf(stderr, "Map of coords from update Coverage:\n");
  //printMap(stderr, mc);
  //fprintf(stderr, "\n------\n");


  map<long, circMapGenomePos>::iterator itp;
  map<int,int>::iterator itmc;
  for(itmc=mc.begin() ; itmc!=mc.end() ; itmc++){
    int posRead = itmc->first;
    int posGenomic = itmc->second;

    //if(chrID=="I" && posGenomic==142619){ fprintf(stderr, "%s\n", readID.c_str()); }

    if(posGenomic != -1){
      char rBase = toupper(readSeq[posRead]);
      if( rBase!='A' && rBase!='T' && rBase!='C' && rBase!='G' )
	continue;
      int iBase = tabCor[rBase];

      //if(chrID=="I" && posGenomic==142619){ fprintf(stderr, "rBase:%c\n", rBase); }	

      itp = (itc->second).find(posGenomic);
      if( itp==(itc->second).end() ){
	simplePosCov spc(rBase);
	spc.add(iBase);
	circMapGenomePos gp;
	gp.set(spc);
	(itc->second).insert( map<long,circMapGenomePos>::value_type(posGenomic,gp) );
	//if(chrID=="I" && posGenomic==142619){ fprintf(stderr, "INSERTING:%c\n", rBase); }	
      } else {
	//if(chrID=="I" && posGenomic==142619){ fprintf(stderr, "UPDATING WITH:%c (current is %c)\n", rBase, itp->second.spc_.getRefBase()); }	
	(itp->second).spc_.add(iBase);
      }
    }
  }
}


void circMapGenome::getRefAndAlterCounts_cov(const std::string &chrID, const int gPos, int *nbRef, int *nbAlter, const char *tabBases){
  map<string, map<long,circMapGenomePos> >::iterator itc = this->mg_.find(chrID);
  this->getRefAndAlterCounts_cov(itc, gPos, nbRef, nbAlter, tabBases);
}


void circMapGenome::getRefAndAlterCounts_cov(map<string, map<long,circMapGenomePos> >::iterator itc, const int gPos, int *nbRef, int *nbAlter, const char *tabBases){
  if( itc==this->mg_.end() ){
    (*nbRef) = 0;
    (*nbAlter) = 0;
  } else {
    map<long,circMapGenomePos>::iterator itp = (itc->second).find(gPos);
    if( itp==(itc->second).end() ){
      (*nbRef) = 0;
      (*nbAlter) = 0;      
    } else {
      (itp->second).spc_.getRefAndAlterCounts(nbRef, nbAlter, tabBases);
    }
  }
}


bool circMapGenome::isSuspect_cov(const std::string &chrID, const int gPos, const int minCov, const double maxFrac, const char *tabBases, int *nbRef, int *nbTot){
  map<string, map<long,circMapGenomePos> >::iterator itc = this->mg_.find(chrID);
  return this->isSuspect_cov(itc, gPos, minCov, maxFrac, tabBases, nbRef, nbTot);
}


bool circMapGenome::isSuspect_cov(map<string, map<long,circMapGenomePos> >::iterator itc, const int gPos, const int minCov, const double maxFrac, const char *tabBases, int *nbRef, int *nbTot){
  if( itc==this->mg_.end() ){ return true; }
  map<long,circMapGenomePos>::iterator itp = (itc->second).find(gPos);
  if( itp==(itc->second).end() ){ return true; }
  return (itp->second).spc_.isSuspect(minCov, maxFrac, tabBases, nbRef, nbTot);
}


bool circMapGenome::updateAllRefBase(map<string,string> &mg, mapOfStrand &mos, char *tab_cpt){
  map<string, map<long,circMapGenomePos> >::iterator itc;
  for(itc=this->mg_.begin() ; itc!=this->mg_.end() ; itc++){
    string chrID = itc->first;
    map<long,circMapGenomePos>::iterator itp;
    for(itp=(itc->second).begin() ; itp!=(itc->second).end() ; itp++){
      long gPos = itp->first;
      itp->second.initRefBase(mg, chrID, gPos, mos, tab_cpt);
    }
  }
  return true;
}

bool circMapGenome::updateAllRefBase(map<string,string> &mg){
  map<string, map<long,circMapGenomePos> >::iterator itc;
  for(itc=this->mg_.begin() ; itc!=this->mg_.end() ; itc++){
    string chrID = itc->first;
    map<long,circMapGenomePos>::iterator itp;
    for(itp=(itc->second).begin() ; itp!=(itc->second).end() ; itp++){
      long gPos = itp->first;
      char refBase = mg[chrID].at(gPos);
      itp->second.setRefBase(refBase);
      //itp->second.initRefBase(mg, chrID, gPos, mos, tab_cpt);
    }
  }
  return true;
}


// END FUNCTIONS THAT DEAL WITH COVERGE

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  *****************************************************************************************************************************************//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// BEGIN FUNCTIONS THAT DEAL WITH MAP OF COVERAGE


void circMapGenome::loadMapIndels(char *fname, bool add){

  string line, chrID, stmp;
  int cov, nbIns, nbDel;
  long gPos;
  string lastChrID = "";

  map<string, map<long,circMapGenomePos> >::iterator itc;
  map<long, circMapGenomePos>::iterator itp;

  ifstream f(fname);
  while( getline(f, line) ){

    istringstream iss(line);
    getline(iss, chrID, '\t');

    if( chrID!=lastChrID ){
      itc = this->mg_.find(chrID);
      if( itc==this->mg_.end() ){
	map<long,circMapGenomePos> mtmp;
	itc = (this->mg_.insert( map<string, map<long,circMapGenomePos> >::value_type(chrID,mtmp) )).first;
      }
      lastChrID = chrID;
    }

    getline(iss, stmp, '\t');
    gPos = atol(stmp.c_str());

    getline(iss, stmp, '\t');
    cov = atoi(stmp.c_str());
    getline(iss, stmp, '\t');
    nbIns = atoi(stmp.c_str());
    getline(iss, stmp, '\t');
    nbDel = atoi(stmp.c_str());

    itp = (itc->second).find(gPos);
    if( itp == (itc->second).end() ){
      C_indelPos ip(cov, nbIns, nbDel);
      (itc->second).insert( map<long,circMapGenomePos>::value_type(gPos, ip) );
    } else {
      if(add==false){
	(itp->second).ip_.set(cov, nbIns, nbDel);
      } else {
	(itp->second).ip_.add(cov, nbIns, nbDel);
      }
      //(itp->second).ip_.add(ip);
    }
  }
}


void circMapGenome::printMapIndels(FILE *nf){
  map<string, map<long,circMapGenomePos> >::iterator itc;
  for(itc = this->mg_.begin() ; itc != this->mg_.end() ; itc++){
    string chrID = itc->first;
    map<long,circMapGenomePos>::iterator itp;
    for(itp=(itc->second).begin() ; itp!=(itc->second).end() ; itp++ ){
      long gPos = itp->first;
      C_indelPos ip = (itp->second).ip_;
      fprintf(nf, "%s\t%d\t", chrID.c_str(), gPos);
      ip.print(nf);
      fprintf(nf, "\n");      
    }
  }
}

int circMapGenome::updateMapIndels(char *samFile){

  samRead sr;
  ifstream fileSam(samFile);

  string previousChrID = "";
  map<string, map<long,circMapGenomePos> >::iterator itc;

  while( sr.initWithNext(fileSam) ){
    string chrID = sr.getRefName();
    itc = this->mg_.find(chrID);

    if( itc == this->mg_.end() ){
      map<long,circMapGenomePos> mtmp;
      itc = (this->mg_.insert( map<string, map<long,circMapGenomePos> >::value_type(chrID, mtmp) )).first;
    }

    map<int,int> mc;
    sr.getMapOfCoords(mc);
    map<int, pair<int,int> > mIns;
    map<int, pair<int,int> > mDel;
    int nbIndel = sr.getIndels(mIns, mDel, mc);
    this->updateCovIndels(itc, mc);
    if(nbIndel>0){
      this->updateIndels(itc, mIns, mDel);
    }
  }
  return 0;
}

int circMapGenome::updateMapIndels(samRead &sr){

  string chrID = sr.getRefName();
  map<string, map<long,circMapGenomePos> >::iterator itc = this->mg_.find(chrID);

  if( itc == this->mg_.end() ){
    map<long,circMapGenomePos> mtmp;
    itc = (this->mg_.insert( map<string, map<long,circMapGenomePos> >::value_type(chrID, mtmp) )).first;
  }

  this->updateMapIndels(sr, itc);
}


int circMapGenome::updateMapIndels(samRead &sr, map<string, map<long,circMapGenomePos> >::iterator itc){

  //string qID = sr.getQueryName();
  //fprintf(stderr, "%s\n", qID.c_str());
  //fflush(stderr);

  map<int,int> mc;
  sr.getMapOfCoords(mc);
  //fprintf(stderr, "map of coords from updateMapIndels:\n");
  //printMap(stderr, mc);
  //fprintf(stderr, "\n--------\n");
  map<int, pair<int,int> > mIns;
  map<int, pair<int,int> > mDel;
  int nbIndel = sr.getIndels(mIns, mDel, mc);

  this->updateCovIndels(itc, mc);
  if(nbIndel>0){
    this->updateIndels(itc, mIns, mDel);
  }

  return 0;
}




int circMapGenome::updateIndels(string &chrID, map<int,pair<int, int> > &mIns, map<int, pair<int,int> > &mDel){
  map<string, map<long,circMapGenomePos> >::iterator itc = this->mg_.find(chrID);
  if( itc == this->mg_.end() ){
    map<long,circMapGenomePos> mtmp;
    itc = (this->mg_.insert( map<string, map<long,circMapGenomePos> >::value_type(chrID,mtmp) )).first;
  }
  return this->updateIndels(itc, mIns, mDel);
}

int circMapGenome::updateIndels(map<string, map<long, circMapGenomePos> >::iterator itc, map<int,pair<int, int> > &mIns, map<int, pair<int,int> > &mDel){

  map<long, circMapGenomePos>::iterator itp;

  map<int, pair<int,int> >::iterator iti;
  for(iti=mIns.begin() ; iti!=mIns.end() ; iti++){
    pair<int,int> pp = iti->second;
    int gPos = pp.first;

    itp = (itc->second).find(gPos);
    if( itp==(itc->second).end() ){
      C_indelPos cp(0,1,0);
      circMapGenomePos gp(cp);
      itp = ( (itc->second).insert( map<long,circMapGenomePos>::value_type(gPos,gp) ) ).first;
    } else {
      (itp->second).ip_.addIns(1);
    }
  }

  map<int, pair<int,int> >::iterator itd;
  for(itd=mDel.begin() ; itd!=mDel.end() ; itd++){
    pair<int,int> pp = itd->second;
    int gPos = pp.first;

    itp = (itc->second).find(gPos);
    if( itp==(itc->second).end() ){
      C_indelPos cp(0,0,1);
      circMapGenomePos gp(cp);
      (itc->second).insert( map<long,circMapGenomePos>::value_type(gPos,gp) );
    } else {
      (itp->second).ip_.addDel(1);
    }
  }

  return 0;
}


int circMapGenome::updateCovIndels(string &chrID, map<int,int> &mc){
  map<string, map<long,circMapGenomePos> >::iterator itc = this->mg_.find(chrID);
  if( itc==this->mg_.end() ){
    map<long,circMapGenomePos> mtmp;
    itc = ( this->mg_.insert( map<string, map<long,circMapGenomePos> >::value_type(chrID,mtmp) ) ).first;
  }
  return this->updateCovIndels(itc, mc);
}


int circMapGenome::updateCovIndels(map<string, map<long,circMapGenomePos> >::iterator itc, map<int,int> &mc){

  if( itc==this->mg_.end() ){ // This should never happen
    fprintf(stderr, "SHOULD NEVER HAPPEN !!!\n");
    fflush(stderr);
    return 1;
  }

  map<long,circMapGenomePos>::iterator itp;

  string chrID = itc->first;

  int firstPos = mc.begin()->second;
  map<int,int>::iterator it;
  for(it=mc.begin() ; it!=mc.end() ; it++){
    int gPos = it->second;
    if(gPos != -1){
      itp = (itc->second).find(gPos);
      if(itp==(itc->second).end()){
	C_indelPos cp(1,0,0);
	circMapGenomePos gp(cp);
	(itc->second).insert( map<long,circMapGenomePos>::value_type(gPos,gp) );
	//if(chrID=="I" && gPos==142619){ fprintf(stderr, "INSERTING:%c (%d)\n", 'N', firstPos); }
      } else {
	itp->second.ip_.addCov(1);
      }
    }
  }
  return 0;
}


bool circMapGenome::isSuspect_indels(std::string chrID, const long gPos, const long minCoverage, const double maxFracAlter){
  map<string, map<long,circMapGenomePos> >::iterator itc = this->mg_.find(chrID);
  if( itc == this->mg_.end() ){ return true; }
  return this->isSuspect_indels(itc, gPos, minCoverage, maxFracAlter);
}


bool circMapGenome::isSuspect_indels(map<string, map<long,circMapGenomePos> >::iterator itc, const long gPos, const long minCoverage, const double maxFracAlter){

  if( itc == this->mg_.end() ){ return true; }

  map<long,circMapGenomePos>::iterator itp = (itc->second).find(gPos);
  if( itp == (itc->second).end() ){ return true; }

  C_indelPos cp = itp->second.ip_;
  return cp.isSuspect(minCoverage, maxFracAlter);
}

// END FUNCTIONS THAT DEAL WITH MAP OF INDELS



// BEGIN FUNCTIONS THAT DEAL WITH MAP OF OBSERVATIONS

void circMapGenome::loadObservations(char *fname, int nbLinesOfHeader, bool erasePreviousMap){

  string line, chrID, sPos, sObs;
  long gPos, nbObs;

  map<string, map<long,circMapGenomePos> >::iterator itc;
  map<long,circMapGenomePos>::iterator itp;

  // First, I set all observations to zero.
  if(erasePreviousMap==true){
    for(itc=this->mg_.begin() ; itc!=this->mg_.end() ; itc++){
      for(itp=(itc->second).begin() ; itp!=(itc->second).end() ; itp++){
	(itp->second).nbObs_ = 0;
      }
    }
  }


  string previousChrID = "";
  ifstream fic(fname);
  for(int i=0 ; i<nbLinesOfHeader ; i++){
    getline(fic,line);
  }
  while( getline(fic, line) ){
    istringstream iss(line);

    getline(iss, chrID, '\t');
    getline(iss, sPos, '\t');
    getline(iss, sObs, '\t');

    gPos = atol(sPos.c_str());
    nbObs = atol(sObs.c_str());

    if( chrID!=previousChrID ){
      itc = this->mg_.find(chrID);
      if( itc==this->mg_.end() ){
	map<long,circMapGenomePos> mtmp;
	itc = ( this->mg_.insert( map<string, map<long,circMapGenomePos> >::value_type(chrID,mtmp) ) ).first;
      }
      previousChrID = chrID;
    }

    itp = (itc->second).find(gPos);
    if( itp==(itc->second).end() ){
      circMapGenomePos gp(nbObs);
      (itc->second).insert( map<long,circMapGenomePos>::value_type(gPos,gp) );
    } else {
      if(erasePreviousMap==true){
	itp->second.nbObs_ = nbObs;
      } else {
	itp->second.nbObs_ += nbObs;
      }
    }
  }  
  fic.close();
}

void circMapGenome::addObservations(char *fname, int nbLinesHeader){
  this->loadObservations(fname, nbLinesHeader, false);
}


void circMapGenome::printObservations(FILE *nf, bool printZeros){
  map<string, map<long,circMapGenomePos> >::iterator itc;
  for(itc=this->mg_.begin() ; itc!=this->mg_.end() ; itc++){
    string chrID = itc->first;
    map<long,circMapGenomePos>::iterator itp;
    for(itp=(itc->second).begin() ; itp!=(itc->second).end() ; itp++){
      long gPos = itp->first;
      long nbObs = itp->second.nbObs_;
      if( printZeros==true || nbObs>0 ){
	fprintf(nf, "%s\t%d\t%ld\n", chrID.c_str(), gPos, nbObs);
      }
    }
  }
}


void circMapGenome::printObservationsWithBase(FILE *nf, map<string,string> &mg, mapOfStrand &mos, bool printStrand=false){

  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);

  map<string, map<long,circMapGenomePos> >::iterator itc;
  for(itc=this->mg_.begin() ; itc!=this->mg_.end() ; itc++){
    string chrID = itc->first;
    map<long,circMapGenomePos>::iterator itp;
    for(itp=(itc->second).begin() ; itp!=(itc->second).end() ; itp++){
      long gPos = itp->first;
      long nbObs = itp->second.nbObs_;

      char genomicBase = mg[chrID].at(gPos);
      char cStrand = mos.getStrand(chrID, gPos);
      if(cStrand!='+' && cStrand!='-'){
	genomicBase = '?';
      } else {
	if(cStrand=='-'){
	  genomicBase = tab_cpt[genomicBase];
	}
      }

      fprintf(nf, "%s\t%d\t%ld\t%c", chrID.c_str(), gPos, nbObs, genomicBase);
      if(printStrand==true){
	fprintf(nf, "\t%c", cStrand);
      }
      fprintf(nf, "\n");
    }
  }
}




int circMapGenome::addObservation(std::string chrID, long gPos){
  map<string, map<long,circMapGenomePos> >::iterator itc = this->mg_.find(chrID);
  if( itc == this->mg_.end() ){
    map<long,circMapGenomePos> mtmp;
    itc = ( this->mg_.insert( map<string, map<long,circMapGenomePos> >::value_type(chrID,mtmp) ) ).first;
  }
  return(this->addObservation(itc, gPos));
  //  (this->mObs_[chrID].at(pos))++;
}

int circMapGenome::addObservation(map<string, map<long,circMapGenomePos> >::iterator itc, long gPos){
  if(itc==this->mg_.end()){ return -1; }

  map<long,circMapGenomePos>::iterator itp = (itc->second).find(gPos);
  if( itp==(itc->second).end() ){
    circMapGenomePos gp(1);
    (itc->second).insert( map<long,circMapGenomePos>::value_type(gPos,gp) );
  } else {
    itp->second.nbObs_++;
  }
}


int circMapGenome::addObservation(std::string chrID, vector<long> & vpos){
  map<string, map<long,circMapGenomePos> >::iterator itc = this->mg_.find(chrID);
  if( itc==this->mg_.end() ){
    map<long,circMapGenomePos> mtmp;
    itc = (this->mg_.insert( map<string, map<long,circMapGenomePos> >::value_type(chrID,mtmp) )).first;
  }

  return this->addObservation(itc, vpos);
}

int circMapGenome::addObservation(map< string,map<long,circMapGenomePos> >::iterator itc, vector<long> & vpos){

  if( itc==this->mg_.end() ){ return -1; }

  vector<long>::iterator itv;
  for(itv=vpos.begin() ; itv!=vpos.end() ; itv++){
    long gPos = (*itv);
    map<long,circMapGenomePos>::iterator itp = (itc->second).find(gPos);
    if( itp == (itc->second).end() ){
      circMapGenomePos gp(1);
      (itc->second).insert( map<long,circMapGenomePos>::value_type(gPos,gp) );
    } else {
      itp->second.nbObs_++;
    }
  }

  return 0;
}

/*
void circMapObs::setObs(std::string chrID, int pos, long nb){
  this->mObs_[chrID].at(pos) = nb;
}

long circMapObs::getObs(std::string chrID, int pos){
  return this->mObs_[chrID].at(pos);
}


long circMapObs::getSum(void){

  long sum = 0;

  map<string, vector<long> >::iterator itm;
  for(itm=this->mObs_.begin() ; itm!=this->mObs_.end() ; itm++){
    vector<long>::iterator itv;
    for(itv=(itm->second).begin() ; itv!=(itm->second).end() ; itv++){
      sum += (*itv);
    }
  }
  return sum;
}
*/

long circMapGenomePos::getNbObservations(void){
  return this->nbObs_;
}


void printMap(FILE *nf, map<int,int> &m){
  map<int,int>::iterator it;
  for(it=m.begin() ; it!=m.end() ; it++){
    fprintf(nf, "%d -> %d\n", it->first, it->second);
  }
}



// END FUNCTIONS THAT DEAL WITH MAP OF OBSERVATIONS


#endif
