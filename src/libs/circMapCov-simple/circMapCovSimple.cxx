#ifndef _CIRC_MAP_COV_SIMPLE_CXX_
#define _CIRC_MAP_COV_SIMPLE_CXX_

#include <stdlib.h>

#include "circMapCovSimple.hxx"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                      simpleCov functions                                                //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

simplePosCov::simplePosCov(void){
  this->init('N');
}

simplePosCov::simplePosCov(const char base){
  this->init(base);
}

int simplePosCov::init(const char base){
  this->refBase_ = base;
  for(int i=0 ; i<4 ; i++){
    this->cov_[i] = 0;
  }
  return 0;
}

int simplePosCov::add(const int iBase){
  this->cov_[iBase]++;
  return 0;
}

int simplePosCov::add(const char base, const int *tabCor){
  int iBase = tabCor[base];
  return this->add(iBase);
}


int simplePosCov::getCoverage(const int iBase){
  return this->cov_[iBase];
}

int simplePosCov::getCoverage(const char base, const int *tabCor){
  int iBase = tabCor[base];
  return this->getCoverage(iBase);
}

char simplePosCov::getRefBase(void){
  return this->refBase_;
}

void simplePosCov::setRefBase(char refBase){
  this->refBase_ = refBase;
}

void simplePosCov::getRefAndAlterCounts(int *nbRef, int *nbAlter, const char* tabBases){
  int iRef=0;
  int iAlter=0;
  for(int i=0 ; i<4 ; i++){
    if(tabBases[i]==this->refBase_){
      iRef = this->cov_[i];
    } else {
      iAlter += this->cov_[i];
    }
  }
  (*nbRef) = iRef;
  (*nbAlter) = iAlter;
}


bool simplePosCov::isSuspect(const int minCov,const double maxFrac, const char *tabBases, int *nbRef, int *nbTot){
  int local_nbRef = 0;
  int local_nbTot = 0;
  int local_nbAlter = 0;
  for(int i=0 ; i<4 ; i++){
    if( tabBases[i] == this->refBase_ ){
      local_nbRef = this->cov_[i];
    } else {
      local_nbAlter += this->cov_[i];
    }
  }

  local_nbTot = local_nbRef + local_nbAlter;
  double frac = (double)((double)local_nbAlter/(double)local_nbTot);

  (*nbRef) = local_nbRef;
  (*nbTot) = local_nbTot;

  if( local_nbTot < minCov || frac > maxFrac )
    return true;

  return false;

}


void simplePosCov::print(FILE *nf) const{
  fprintf(nf, "%c", this->refBase_);
  for(int i=0 ; i<4 ; i++){
    fprintf(nf, "\t%d", this->cov_[i]);
  }
}


int simplePosCov::set(const int iBase, const int nb){
  this->cov_[iBase] = nb;
}

int simplePosCov::set(const char base, const int nb, const int*tabCor){
  int iBase = tabCor[base];
  return this->set(iBase, nb);
}


void simplePosCov::add(simplePosCov &spc){
  for(int i=0 ; i<4 ; i++){
    this->cov_[i] += spc.getCoverage(i);
  }
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                  circMapCovSimple functions                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

circMapCovSimple::circMapCovSimple(void){
  ;
}

circMapCovSimple::circMapCovSimple(char *fileGenome){
  map<string,string> mg;
  loadFasta(mg, fileGenome);
  this->init(mg);
}

circMapCovSimple::circMapCovSimple(char *fileGenome, char *fileSam, const int *tabCor, const int minScore){
  map<string,string> mg;
  loadFasta(mg, fileGenome);
  this->init(mg);
  this->update(fileSam, tabCor, minScore);
}


circMapCovSimple::circMapCovSimple(const std::map<std::string,std::string> &mg){
  this->init(mg);
}


void circMapCovSimple::load(char *fname){

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

  map<string, map<long, simplePosCov> >::iterator it;
  string prevChrID = "";

  while( getline(fic, line) ){
    istringstream iss(line);
    getline(iss, chrID, '\t');
    getline(iss, stmp, '\t');
    long gPos = atol(stmp.c_str());

    getline(iss, stmp, '\t');
    char refBase = stmp[0];
    simplePosCov spc(refBase);

    if(chrID!=prevChrID){
      it = this->mc_.find(chrID);
      if(it==this->mc_.end()){
	map<long, simplePosCov> m_tmp;
	this->mc_[chrID] = m_tmp;
	it = this->mc_.find(chrID);
      }
      prevChrID = chrID;
    }

    for(i=0 ; i<4 ; i++){
      getline(iss, stmp, '\t');
      int nb = atoi(stmp.c_str());
      spc.set(i,nb);
      //(it->second).at(gPos).set(i, nb);
      //this->mc_[chrID].at(gPos).set(i, nb);
    }

    (it->second).insert( map<long,simplePosCov>::value_type(gPos,spc) );

  }

}


void circMapCovSimple::load(char *fname, string chrOnly){

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

  map<string, map<long,simplePosCov> >::iterator it;
  string prevChrID = "";

  while( getline(fic, line) ){
    istringstream iss(line);
    getline(iss, chrID, '\t');

    if(chrID!=chrOnly)
      continue;

    getline(iss, stmp, '\t');
    long gPos = atol(stmp.c_str());

    getline(iss, stmp, '\t');
    char refBase = stmp[0];
    simplePosCov spc(refBase);

    if(chrID!=prevChrID){
      it = this->mc_.find(chrID);
      if(it==this->mc_.end()){
	map<long, simplePosCov> m_tmp;
	this->mc_[chrID] = m_tmp;
	it = this->mc_.find(chrID);
      }
      prevChrID = chrID;
    }

    for(i=0 ; i<4 ; i++){
      getline(iss, stmp, '\t');
      int nb = atoi(stmp.c_str());
      spc.set(i, nb);
      //(it->second).at(gPos).set(i, nb);
      //this->mc_[chrID].at(gPos).set(i, nb);
    }

    (it->second).insert( map<long,simplePosCov>::value_type(gPos,spc) );

  }

}



int circMapCovSimple::init(const map<string,string> &mg){

  map<string,string>::const_iterator itmg;

  for(itmg=mg.begin() ; itmg!=mg.end() ; itmg++){
    string chrID = itmg->first;
    string seq = itmg->second;
    int lg = seq.size();

    map<long, simplePosCov> m_tmp;
    this->mc_[chrID] = m_tmp;
  }
  return 0;
}


int circMapCovSimple::update(const char *fileSam, const int *tabCor, const int minScore){

  samRead sr;
  ifstream ficSam(fileSam);

  while( sr.readNext(ficSam) ){
    this->processSamRead(sr, tabCor, minScore);
  }

}


int circMapCovSimple::processSamRead(const samRead &sr, const int* tabCor, const int minScore){
  map<int,int> mc;
  sr.getGenomicCoords(mc);

  string chrID = sr.getRefName();
  string readSeq = sr.getSeq();

  map<string, map<long,simplePosCov> >::iterator itm = this->mc_.find(chrID);
  map<long, simplePosCov>::iterator itms;

  map<int,int>::iterator itmc;
  for(itmc=mc.begin() ; itmc!=mc.end() ; itmc++){
    int posRead = itmc->first;
    int posGenomic = itmc->second;

    if(posGenomic != -1){
      char rBase = toupper(readSeq[posRead]);
      if( rBase!='A' && rBase!='T' && rBase!='C' && rBase!='G' )
	continue;
      int iBase = tabCor[rBase];

      itms = (itm->second).find(posGenomic);
      if( itms==(itm->second).end() ){
	simplePosCov spc(rBase);
	spc.add(iBase);
	(itm->second).insert( map<long,simplePosCov>::value_type(posGenomic,spc) );
      } else {
	(itms->second).add(iBase);
      }

    }
  }
}


void circMapCovSimple::getRefAndAlterCounts(const std::string &chrID, const int gPos, int *nbRef, int *nbAlter, const char *tabBases){
  map<string, map<long,simplePosCov> >::iterator it = this->mc_.find(chrID);
  this->getRefAndAlterCounts(it, gPos, nbRef, nbAlter, tabBases);
}


void circMapCovSimple::getRefAndAlterCounts(map<string, map<long,simplePosCov> >::iterator it, const int gPos, int *nbRef, int *nbAlter, const char *tabBases){

  if( it==this->mc_.end() ){
    (*nbRef) = 0;
    (*nbAlter) = 0;
  } else {
    map<long,simplePosCov>::iterator itm = (it->second).find(gPos);
    if( itm==(it->second).end() ){
      (*nbRef) = 0;
      (*nbAlter) = 0;      
    } else {
      (itm->second).getRefAndAlterCounts(nbRef, nbAlter, tabBases);
    }
  }
}



bool circMapCovSimple::isSuspect(const std::string &chrID, const int gPos, const int minCov, const double maxFrac, const char *tabBases, int *nbRef, int *nbTot){

  map<string, map<long,simplePosCov> >::iterator it = this->mc_.find(chrID);
  return this->isSuspect(it, gPos, minCov, maxFrac, tabBases, nbRef, nbTot);

}


bool circMapCovSimple::isSuspect(map<string, map<long,simplePosCov> >::iterator it, const int gPos, const int minCov, const double maxFrac, const char *tabBases, int *nbRef, int *nbTot){

  if( it==this->mc_.end() ){
    return true;
  }

  map<long,simplePosCov>::iterator itm = (it->second).find(gPos);

  if( itm==(it->second).end() ){
    return true;
  }

  return (itm->second).isSuspect(minCov, maxFrac, tabBases, nbRef, nbTot);

}


void circMapCovSimple::print(FILE *nf, const char *tabBases) const{

  int i;

  fprintf(nf, "%s\t%s\t%s", "scafID", "pos", "refBase");
  for(i=0 ; i<4 ; i++){
    fprintf(nf, "\t%c", tabBases[i]);
  }
  fprintf(nf, "\n");

  map<string, map<long,simplePosCov> >::const_iterator it;
  for(it=this->mc_.begin() ; it!=this->mc_.end() ; it++){
    string chrID = it->first;

    map<long, simplePosCov> mm = (it->second);
    map<long, simplePosCov>::iterator itm;

    for(
	itm = mm.begin() ;
	itm != mm.end() ;
	itm++
	){
      long iPos = itm->first;
      simplePosCov sc = itm->second;
      fprintf(nf, "%s\t%d\t", chrID.c_str(), iPos);
      sc.print(nf);
      fprintf(nf, "\n");
    }
  }

}

#endif
