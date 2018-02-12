#ifndef _SAM_READ_CXX_
#define _SAM_READ_CXX_


#include <stdlib.h>

#include "samRead.hxx"


using namespace std;


void samRead_vectoriseCigar(vector< pair<char, int> > &v, const string &s);



samRead::~samRead(void){
  this->mOpt_.clear();
}

samRead::samRead(void){
  ;
}

samRead::samRead(string line){
  this->init(line);
}



int samRead::init(string line){

  string stmp;
  this->mOpt_.clear();

  istringstream iss(line);

  getline(iss, this->qName_, '\t');

  getline(iss, stmp, '\t');
  this->flag_ = atoi(stmp.c_str());

  getline(iss, this->rName_, '\t');

  getline(iss, stmp, '\t');
  this->pos_ = atoi(stmp.c_str()) - 1;

  getline(iss, stmp, '\t');
  this->mapq_ = atoi(stmp.c_str());

  getline(iss, this->cigar_, '\t');
  getline(iss, this->rNext_, '\t');

  getline(iss, stmp, '\t');
  if( this->rNext_ == "*" ){
    this->pNext_ = -1;
  } else {
    this->pNext_ = atoi(stmp.c_str()) - 1;
  }

  getline(iss, stmp, '\t');
  if( this->rNext_ == "*" ){
    this->tLen_ = -1;
  } else {
    this->tLen_ = atoi(stmp.c_str());
  }


  getline(iss, this->seq_, '\t');
  getline(iss, this->qual_, '\t');


  /**/
  string sKey, sType, sVal;
  while( getline(iss, stmp, '\t') ){
    istringstream iss2(stmp);
    getline(iss2, sKey, ':');
    getline(iss2, sType, ':');
    getline(iss2, sVal, ':');
    this->mOpt_.insert(map<string,string>::value_type(sKey, sVal));
  }
  /**/

}

bool samRead::initWithNext(ifstream &f){
  string line;

  while( getline(f, line) ){

    if( line.size()>0 && line[0]!='@' ){
      this->init(line);
      return true;
    }

  }

  return false;
}


int samRead::readNext(std::ifstream &fic){
  string line;

  if( getline(fic, line) ){

    if( line.size()>0 && line[0]=='@' ) // I'm looking at a header line
      return -1;

    this->init(line);
    return 1;
  }

  return 0;
}


int samRead::readNext(std::ifstream &fic, string &line){
  if( getline(fic, line) ){

    if( line.size()>0 && line[0]=='@' )// I'm looking at a header line
      return -1;

    this->init(line);
    return 1;
  }

  return 0;
}

int samRead::getNbIndels(void){
  int nbIndel = 0;
  vector< pair<int,char> > vsc;
  this->getSplitCigar(vsc);
  int nb = vsc.size();
  for(int i=0 ; i<nb ; i++){
    char cType = vsc[i].second;
    if(cType=='I' || cType=='D'){
      nbIndel++;
    }
  }
  return nbIndel;
}


bool samRead::getGenomicCoords(map<int,int> &m) const{

  vector< pair<char, int> > vCigar;
  samRead_vectoriseCigar(vCigar, this->cigar_);

  int posGenomic = this->pos_;
  int posRead = 0;
  int ii;

  vector< pair<char,int> >::iterator itv;
  for(itv=vCigar.begin() ; itv!=vCigar.end() ; itv++){
    char code = (*itv).first;
    int nb = (*itv).second;

    // I do not know how to handle the (rare) reads with exotic codes in their CIGAR string
    if(code!='M' && code!='N' && code!='I' && code!='D'){
      m.clear();
      return false;
    }

    if(code=='M'){
      for(ii=0 ; ii<nb ; ii++, posRead++){
	m[posRead] = posGenomic + ii;
      }
      posGenomic += nb;
    }

    if(code=='N' || code=='D'){
      posGenomic += (nb); // <- XX this +1 might be extra and create a shift...
    }

    if(code=='I'){
      for(ii=0 ; ii<nb ; ii++, posRead++){
	m[posRead] = -1;
      }
    }
  }

  return true;
}

// !!! To DO !!! -- STILL NEEDS TO BE FINISHED
//int samRead::getIndelsCoords(map<int,int> &mIns, map<int,int> &mDel) const{

  // HERE: INTENTIONAL COMPILATION ERROR.
  /*
  if( this->cigar_.find('D')==string::npos && this->cigar_.find('I')==string::npos ){
    return 0;
  }

  vector< pair<char, int> > vCigar;
  samRead_vectoriseCigar(vCigar, this->cigar_);

  int posGenomic = this->pos_;
  int posRead = 0;
  int ii;

  vector< pair<char,int> >::iterator itv;
  for(itv=vCigar.begin() ; itv!=vCigar.end() ; itv++){
    char code = (*itv).first;
    int nb = (*itv).second;

    // I do not know how to handle the (rare) reads with exotic codes in their CIGAR string
    if(code!='M' && code!='N' && code!='I' && code!='D'){
      return -1;
    }

    if(code=='M' || code=='N'){
      posGenomic += nb;
    }

    if(code=='D'){
      
      posGenomic += nb;
    }

    if(code=='I'){
      for(ii=0 ; ii<nb ; ii++, posRead++){
	m[posRead] = -1;
      }
    }
  }
*/
//  return 1;
//}


string samRead::getGenomicSeq(map<int,int> &mc, map<string, string> &mg) const{

  map<int,int>::iterator itmc;

  string chrID = this->rName_;
  string seqGenomic = "";

  for(itmc=mc.begin() ; itmc!=mc.end() ; itmc++){
    int posRead = itmc->first;
    int posGenomic = itmc->second;
    if( posGenomic == -1 ){
      seqGenomic += "-";
    } else {
      seqGenomic += mg[chrID].at(posGenomic);
    }
  }

  return seqGenomic;
}


bool samRead::isPaired(void){ return (this->flag_ & 0x1); }
bool samRead::readMappedInProperPair(void){ return (this->flag_ & 0x2); }
bool samRead::readUnmapped(void){ return (this->flag_ & 0x4); }
bool samRead::mateUnmapped(void){ return (this->flag_ & 0x8); }
bool samRead::readReverseStrand(void){ return (this->flag_ & 0x10); }
bool samRead::mateReverseStrand(void){ return (this->flag_ & 0x20); }
bool samRead::firstInPair(void){ return (this->flag_ & 0x40); }
bool samRead::secondInPair(void) { return (this->flag_ & 0x80); }
bool samRead::notPrimaryAlignment(void){ return (this->flag_ & 0x100); }
bool samRead::readFailsQualityChecks(void){ return (this->flag_ & 0x200); }
bool samRead::readIsPCRduplicate(void){ return (this->flag_ & 0x400); }
bool samRead::supplementaryAlignment(void){ return (this->flag_ & 0x800); }


string samRead::getQueryName(void) const{
  return this->qName_;
}

string samRead::getRefName(void) const{
  return this->rName_;
}

string samRead::getCigar(void) const{
  return this->cigar_;
}

string samRead::getSeq(void) const{
  return this->seq_;
}

int samRead::getSeqLg(void) const{
  return this->seq_.size();
}

int samRead::getFlag(void) const{
  return this->flag_;
}

int samRead::getNbMM(void) {
  string sKey = "XM";
  string sVal = this->mOpt_[sKey];
  int nbMM = atoi( sVal.c_str() );
  return nbMM;
}


int samRead::getEditDistanceFromFlag(void){
  string sKey = "NM";
  string sVal = this->mOpt_[sKey];
  int dist = atoi(sVal.c_str());
  return dist;
}

int samRead::getEditDistance(void){
  int nbMM, nbIndels, nbIns, nbDel, totalSizeIndels, totalSizeIns, totalSizeDel, nbN, nbOthers;
  int dist = this->getEditDistance(&nbMM, &nbIndels, &nbIns, &nbDel, &totalSizeIndels, &totalSizeIns, &totalSizeDel, &nbN, &nbOthers);
  return dist;
}


int samRead::getEditDistance(int *nbMM, int *nbIndels, int *nbIns, int *nbDel, int *totalSizeIndels, int *totalSizeIns, int *totalSizeDel, int *nbN, int *nbOthers){

  (*nbMM) = this->getNbMM();
  (*nbIndels) = 0;

  (*nbIns) = 0;
  (*nbDel) = 0;
  (*totalSizeIndels) = 0;
  (*totalSizeIns) = 0;
  (*totalSizeDel) = 0;

  (*nbN) = 0;
  (*nbOthers) = 0;

  int nbM = 0;

  vector< pair<int, char> > vsc;
  this->getSplitCigar(vsc);
  int splitSize = vsc.size();

  for(int i=0 ; i<splitSize ; i++){
    char cType = vsc[i].second;
    int cSize = vsc[i].first;
    switch(cType){
    case 'M':
      nbM++;
      break;
    case 'I':
      (*nbIndels)++;
      (*nbIns)++;
      (*totalSizeIndels) += cSize;
      (*totalSizeIns) += cSize;
      break;
    case 'D':
      (*nbIndels)++;
      (*nbDel)++;
      (*totalSizeIndels) += cSize;
      (*totalSizeDel) += cSize;
      break;
    case 'N':
      (*nbN)++;
      break;
    default:
      (*nbOthers)++;
      break;
    }
  }

  int dist = (*nbMM) + (*nbIndels);
  return dist;

}


int samRead::getPos(void){
  return this->pos_;
}

void samRead::setPos(int pos){
  this->pos_ = pos;
}


bool samRead::getSplitCigar(std::vector< std::pair<int, char> > & v) const{

  char cSize[12];
  char car;
  int i, lPos;
  int lg = this->cigar_.size();

  for(i=0, lPos=0 ; i<lg ; i++){
    char c = this->cigar_[i];
    if(c<'0' || c>'9'){
      cSize[lPos]='\0';
      int lSize = atoi(cSize);
      pair<int, char> pp(lSize, c);
      v.push_back(pp);
      lPos = 0;
    } else {
      cSize[lPos] = c;
      lPos++;
    }
  }


}

int samRead::getCigarDepth(void){
  int depth = 0;
  int i;
  int cigarLength = this->cigar_.size();
  for(i=0 ; i<cigarLength ; i++){
    char car = this->cigar_[i];
    if( car<'0' || car>'9' ){ depth++; }
  }
  return depth;
}


/*
int samRead::getIndelPosAndSize(int *pos, int *size){

  vector< pair<int, char> > vsp;
  this->getSplitCigar(vsp);

  int i;
  int nbs = vsp.size();

  // I look only at cases with maximum of 1 indel and this indel has to be between 2 matches.
  if(nbs!=3){
    return -1;
  }

  if( vsp[0].second!='M' || vsp[2].second!='M' ){
    fprintf(stderr, "WARNING, CIGAR STRING WITH INDEL AT EXTREMITY (%s)\n", this->qName_.c_str());
    return -1;
  }

  (*pos) = vsp[0].first;
  (*size) = vsp[1].first;
  int c_indelType = (int)vsp[1].second;
  return c_indelType;

}
*/

int samRead::refine(map<int, int> &mc, string & chrSeq, GFFtranscript &gt, int *offset){

  //this->printAln(stderr, mc, chrSeq);
  //fprintf(stderr, "\n----\n");
  int i;
  int nbMM = this->getNbMM();
  int eDist = this->getEditDistance();

  int nbIndel = 0;
  char indelType = (char)-1;
  int indelPosInRead = -1;
  int indelSize = -1;


  vector< pair<int, char> > vsc;
  this->getSplitCigar(vsc);
  int splitSize = vsc.size();
  int nbM = 0;
  int nbN = 0;
  int nbOthers = 0;
  int nbMbefore = 0;

  for(int i=0 ; i<splitSize ; i++){
    char cType = vsc[i].second;
    int cSize = vsc[i].first;
    switch(cType){
    case 'M':
      nbM++;
      nbMbefore += cSize;
      break;
    case 'I':
      nbIndel++;
      indelType = 'I';
      indelSize = cSize;
      indelPosInRead = nbMbefore - 1;
      break;
    case 'D':
      nbIndel++;
      indelType = 'D';
      indelSize = cSize;
      indelPosInRead = nbMbefore - 1;
      break;
    case 'N':
      nbN++;
      break;
    default:
      nbOthers++;
      break;
    }
  }


  if( nbIndel+nbMM>1 || nbOthers>0 )
    return -1;


  int nbTotEdit = nbIndel + nbMM;
  if( nbTotEdit!= 1 ) // I look only at mapped reads with at most one indels OR one mismatch (--> I reject cases with one indels AND one mismatch)
    return -1;

  //map<int,int> mc;
  //this->getMapOfCoords(mc);

  int editPos;
  if(nbIndel==1){
    editPos = indelPosInRead;
  } else {
    editPos = this->findFirstEdit(mc, chrSeq);
  }

  //fprintf(stderr, "EditPos:%d\n", editPos);

  if(editPos==-1){
    fprintf(stderr, "PROBLEM: findFirstEdit returned -1 (should never happen) - readID for debog: '%s'\n", this->qName_.c_str());
    return -1;
  }
  if(editPos==-2){
    fprintf(stderr, "PROBLEM: findFirstEdit did not find any edit (should never happen) - readID for debog: '%s'\n", this->qName_.c_str());
    return -2;
  }


  int readLg = this->seq_.size();
  map<int,int> mTmp;
  double rPosEdit = (double)editPos/(double)readLg;
  int nbToMove = -1;
  int newDist;

  if( rPosEdit>0.5 ){ // If the mismatch/indel is in the second half of the consensus sequence, then I need to move nucleotides from the end of the consensus to the begining
    if(nbIndel==1){
      nbToMove = (readLg - editPos) - 1;
    } else {
      nbToMove = readLg - editPos;
    }

    //fprintf(stderr, "NbToMove (before):%d\n", nbToMove);
    newDist = this->moveToLeft(mc, mTmp, &nbToMove, gt, chrSeq, indelType, indelPosInRead, indelSize);
    //fprintf(stderr, "NbToMove (after):%d\n", nbToMove);
    (*offset) = (-nbToMove);
  } else {
    nbToMove = editPos + 1;
    if( indelType=='I' ){
      nbToMove += (indelSize);
    }
    //fprintf(stderr, "NbToMove (before):%d\n", nbToMove);
    newDist = this->moveToRight(mc, mTmp, &nbToMove, gt, chrSeq, indelType, indelPosInRead, indelSize);
    //fprintf(stderr, "NbToMove (after):%d\n", nbToMove);
    (*offset) = nbToMove;
  }

  //int newDist = this->getEditDistance(mTmp, chrSeq);
  //(*offset) = 1;
  //fprintf(stderr, "offset:%d\n", (*offset));

  //fprintf(stderr, "After refine:\n");
  //this->printAln(stderr, mc, chrSeq);
  //fprintf(stderr, "\n----\n");

  return newDist;

}


int samRead::moveToRight(map<int,int> &m, map<int,int> &newM, int * nbToMove, GFFtranscript &gt, string &chrSeq, int indelType, int indelPos, int indelSize){

  //this->printAln(stderr, m, chrSeq);

  int i, posRead;
  int seqLg = this->seq_.size();
  int lnb = (*nbToMove);

  int maxPosGenomic = m[seqLg-1];
  if(maxPosGenomic<0){
    fprintf(stderr, "ERROR: maxPosGenomic<0 ('%s')\n", this->qName_.c_str());
    return -1;
  }

  int maxPosTranscript = gt.genomicToSplicedTranscript(maxPosGenomic);

  //////////////////////////////////////////////////////////////////////////////
  // DEBOG
  //for(i=seqLg-1 ; i>=0 ; i--){
  //fprintf(stderr, "%d: %d - %d\n", i, m[i], gt.genomicToSplicedTranscript(m[i]));
  //}


  //if(maxPosTranscript<0){
  //fprintf(stderr, "WARNING: First initialization of maxPosTranscript return -1. Read id for debog: '%s'\n", this->qName_.c_str());
  //}

  int lastPosGenomic = maxPosGenomic;
  int posGenomic = maxPosGenomic;
  int posTranscript = maxPosTranscript;

  bool onPlusStrand = gt.isOnPlusStrand();
  bool mappingNeverInGene = true;

  for(i=0 ; i<lnb ; i++){

    if(posTranscript>=0){
      mappingNeverInGene = false;
      if(onPlusStrand==true){
	posTranscript++;
      } else {
	posTranscript--;
      }    
      mappingNeverInGene = false;
      posGenomic = gt.splicedToGenomic(posTranscript);
      lastPosGenomic = posGenomic;
    } else {
      posGenomic = lastPosGenomic + 1;
      lastPosGenomic = posGenomic;
      posTranscript = gt.genomicToSplicedTranscript(lastPosGenomic);
      if( posTranscript>=0 ){ // Here is the tricky case of when we move from the non-annotated to the annotated part of the mapping.
	// I artificially position the 'posTranscript' variable one nucleotide beyond the end of the transcript so that the next shift will get it to the correct position.
	if(onPlusStrand==true){
	  posTranscript--;
	} else {
	  posTranscript++;
	}
      }
    }

    if( posGenomic<0 ){
      //fprintf(stderr, "ERROR: posGenomic<0 in move to right ('%s')\n", this->qName_.c_str());
      //return -1;
      posGenomic = lastPosGenomic + 1;
      lastPosGenomic = posGenomic;
      posTranscript = gt.genomicToSplicedTranscript(lastPosGenomic);
    }

    newM[i] = posGenomic;
  }

  if(mappingNeverInGene==true){
    //fprintf(stderr, "WARNING, '%s' mapped completely outside any annotated gene!\n", this->qName_.c_str());
    ;
  }

  int newDist = this->getEditDistance(newM, chrSeq);
  return newDist;
}


int samRead::moveToLeft(map<int,int> &m, map<int,int> &newM, int * nbToMove, GFFtranscript &gt, string &chrSeq, int indelType, int indelPos, int indelSize){

  newM = m;
  //myLocal_samRead_printMap(stderr, m);

  int i, posRead;
  int seqLg = this->seq_.size();
  int lnb = (*nbToMove);

  int minPosGenomic = m[0];
  if(minPosGenomic<0){
    fprintf(stderr, "ERROR: minPosGenomic<0 ('%s')\n", this->qName_.c_str());
    return -1;
  }
  int minPosTranscript = gt.genomicToSplicedTranscript(minPosGenomic);
  int test = gt.splicedToGenomic(minPosTranscript);


  //if(minPosTranscript<0){
  //fprintf(stderr, 
  //	    "WARNING: First initialization of minPosTranscript return -1 (posGenomic=%d , transcriptID:'%s'). Read id for debog: '%s'\n", 
  //	    minPosGenomic, 
  //	    gt.id_.c_str(), 
  //	    this->qName_.c_str()
  //	    );
  //}

  int lastPosGenomic = minPosGenomic;
  int posGenomic = minPosGenomic;
  int posTranscript = minPosTranscript;

  bool onPlusStrand = gt.isOnPlusStrand();
  bool mappingNeverInGene = true;


  // In case I'm looking at a deletion, I try first to fill in the deletion as much as possible:
  if( indelType=='D' || indelType=='I' ){
    int nbFill = this->tryToFillDeletionLeft(newM, indelPos, indelSize, gt, chrSeq);
    (*nbToMove) -= nbFill;
    lnb = (*nbToMove);
  }


  // Here, add some code to keep in place all the mapping nucleotides (usefull for example if the deletion is called on the first nucleotide of a homopolymer run, the del should actually be called on the last one...

  for(i=(seqLg-1) ; i>=(seqLg-lnb) ; i--){

    if( posTranscript>=0 ){
      if(onPlusStrand==true){
	posTranscript--;
      } else {
	posTranscript++;
      }
      mappingNeverInGene = false;
      posGenomic = gt.splicedToGenomic(posTranscript);
      lastPosGenomic = posGenomic;
    } else {
      posGenomic = lastPosGenomic - 1;
      lastPosGenomic = posGenomic;
      posTranscript = gt.genomicToSplicedTranscript(lastPosGenomic);
      if( posTranscript>=0 ){ // Here is the tricky case of when we move from the non-annotated to the annotated part of the mapping.
	// I artificially position the 'posTranscript' variable one nucleotide beyond the end of the transcript so that the next shift will get it to the correct position.
	if(onPlusStrand==true){
	  posTranscript++;
	} else {
	  posTranscript--;
	}
      }
    }

    if( posGenomic<0 ){ // This means I'm getting beyond the end of the annotated gene.
      //fprintf(stderr, "ERROR: posGenomic<0 in move to left ('%s')\n", this->qName_.c_str());
      //return -1;
      posGenomic = lastPosGenomic - 1;
      lastPosGenomic = posGenomic;
      posTranscript = gt.genomicToSplicedTranscript(lastPosGenomic);
    }

    newM[i] = posGenomic;
  }

  //myLocal_samRead_printMap(stderr, newM);
  if(mappingNeverInGene==true){
    fprintf(stderr, "WARNING, '%s' mapped completely outside any annotated gene!\n", this->qName_.c_str());
  }

  int newDist = this->getEditDistance(newM, chrSeq);

  if( newDist>0 && (indelType=='D' || indelType=='I') ){
    int nb = 0;
    char bRead, bGenomic;
    bool firstRound = true;
    while(firstRound==true || (nb<indelSize && posGenomic!=-1 && bRead==bGenomic) ){
      firstRound = false;
      int posInRead = (indelPos + nb);
      posTranscript = gt.genomicToSplicedTranscript(newM[posInRead-1]);

      if(onPlusStrand==true){
	posTranscript++;
      } else {
	posTranscript--;
      }

      posGenomic = gt.splicedToGenomic(posTranscript);
      bRead = this->seq_[posInRead];
      bGenomic = '?';
      if(posGenomic>=0){
	bGenomic = chrSeq[posGenomic];
      }

      //fprintf(stderr, "%c:%c (posRead:%d - posTranscript:%d - posGenomic:%d\n", bRead, bGenomic, posInRead, posTranscript, posGenomic);

      if(bRead==bGenomic){
	newM[posInRead] = posGenomic;
	(*nbToMove)--;
      }
      nb++;
    }

  }

  //fprintf(stderr, "newM:\n");
  //this->printAln(stderr, newM, chrSeq);
  //myLocal_samRead_printMap(stderr, newM);

  newDist = this->getEditDistance(newM, chrSeq);
  return newDist;
}



int samRead::getMapOfCoords(map<int, int> &m) const{

  int i, j;

  vector< pair<int, char> > vsc;
  this->getSplitCigar(vsc);
  int splitSize = vsc.size();

  int posGenomic = this->pos_;
  int posInRead = 0;

  for(i=0 ; i<splitSize ; i++){
    char cType = vsc[i].second;
    int cSize = vsc[i].first;

    if( cType!='M' && cType!='N' && cType!='D' && cType!='I' )
      return -1;

    if( cType=='M' ){
      for(j=0 ; j<cSize ; j++){
	m[posInRead] = posGenomic;
	posInRead++;
	posGenomic++;
      }
    }

    if( cType=='D' || cType=='N' ){
      posGenomic += cSize;
    }

    if( cType=='I' ){
      for(j=0 ; j<cSize ; j++){
	m[posInRead] = -1;
	posInRead++;
      }
    }

  }

  return 0;  
}


int samRead::getMapOfCoordsV2(map<int, int> &m){

  int i, j;

  vector< pair<int, char> > vsc;
  this->getSplitCigar(vsc);
  int splitSize = vsc.size();

  int posGenomic = this->pos_;
  int posInRead = 0;

  for(i=0 ; i<splitSize ; i++){
    char cType = vsc[i].second;
    int cSize = vsc[i].first;

    if( cType!='M' && cType!='N' && cType!='D' && cType!='I' )
      return -1;

    if( cType=='M' ){
      for(j=0 ; j<cSize ; j++){
	m[posInRead] = posGenomic;
	posInRead++;
	posGenomic++;
      }
    }

    if( cType=='D' || cType=='N' ){
      posGenomic += cSize;
    }

    if( cType=='I' ){
      for(j=0 ; j<cSize ; j++){
	m[posInRead] = -1;
	posInRead++;
      }
    }

  }

  return 0;  
}





// Same as findFirstMismatch, except that I would also return an indel.
int samRead::findFirstEdit(map<int,int> &m, string &chrSeq){

  int i;

  vector< pair<int, char> > vsc;
  this->getSplitCigar(vsc);
  int splitSize = vsc.size();

  int firstIndelPos = -1;
  int nbMatchBefore = 0;
  bool indelFound = false;
  for(i=0 ; indelFound==false && i<splitSize ; i++){
    char cType = vsc[i].second;
    int cSize = vsc[i].first;

    if(cType=='M'){
      nbMatchBefore += cSize;
    }

    if( cType=='I' || cType=='D' ){
      firstIndelPos = nbMatchBefore;
      indelFound = true;
    }
  }

  int posInRead, posGenomic;
  int chrSize = chrSeq.size();

  map<int,int>::iterator it;
  for(it=m.begin() ; it!=m.end() ; it++){
    posInRead = it->first;
    posGenomic = it->second;

    if(indelFound==true && posInRead>firstIndelPos)
      return firstIndelPos;


    if(posGenomic<0 || posGenomic>=chrSize)
      return -1;

    char bRead = toupper(this->seq_[posInRead]);
    char bGenome = toupper(chrSeq[posGenomic]);
    if(bRead!=bGenome){
      return posInRead;
    }
  }

  return -2;

}


int samRead::findFirstMismatch(map<int,int> &m, string &chrSeq){

  int posInRead, posGenomic;
  map<int, int>::iterator it;

  int nbDiff = 0;
  int chrSize = chrSeq.size();

  for(it=m.begin() ; it!=m.end() ; it++){
    posInRead = it->first;
    posGenomic = it->second;

    if(posGenomic>=0){ // I skip the insertions, they are not considered as 'mismatches'

      if(posGenomic>=chrSize)
	return -1;

      char bRead = toupper(this->seq_[posInRead]);
      char bGenome = toupper(chrSeq[posGenomic]);
      if(bRead!=bGenome){
	return posInRead;
      }

    }
  }

  return -2;
}



int samRead::getEditDistance(map<int,int> & m, string & chrSeq){

  int posInRead, posGenomic;
  map<int, int>::iterator it;

  int nbDiff = 0;
  int chrSize = chrSeq.size();

  for(it=m.begin() ; it!=m.end() ; it++){
    posInRead = it->first;
    posGenomic = it->second;

    if(posGenomic<0 || posGenomic>=chrSize)
      return -1;

    char bRead = toupper(this->seq_[posInRead]);
    char bGenome = toupper(chrSeq[posGenomic]);
    if(bRead!=bGenome){
      nbDiff++;
    }
  }
  return nbDiff;
}


int samRead::getIndels(std::map<int, std::pair<int,int> > &mIns, std::map<int, std::pair<int,int> > &mDel){
  map<int,int> mc;
  this->getMapOfCoords(mc);
  return( this->getIndels(mIns, mDel, mc) );
}


// structure of map<int, pair<int,int> > :
// map< positionIndelInRead , pair< positionInGenome , indelSize > >
// For deletions, the position in the read is the position upstream of the deletion.
// For insertions, the position in the genome is the position upstream of the insertion.
int samRead::getIndels(std::map<int, std::pair<int,int> > &mIns, std::map<int, std::pair<int,int> > &mDel, std::map<int,int> &mc){

  int i = 0;
  int nbIndel = 0;

  int posGenome = this->pos_;
  int posRead = 0;

  mIns.clear();
  mDel.clear();

  vector< pair<int, char> > vs;
  this->getSplitCigar(vs);
  int nb = vs.size();

  if(nb==1){
    return 0;
  }

  for(i=0 ; i<nb ; i++){
    int size = vs[i].first;
    char cType = vs[i].second;    
    switch(cType){
    case 'M':
    case 'N':
      posRead += size;
      posGenome += size;
      break;
    case 'D': {
      pair<int,int> ppd(posGenome-1, size);
      mDel[posRead-1] = ppd;
      posGenome += size;
      nbIndel++;
      break;
    }
    case 'I': {
      pair<int,int> ppi(posGenome, size);
      mIns[posRead] = ppi;
      posRead += size;
      nbIndel++;
      break;
    }
    default:
      break;
    }
  }

  return nbIndel;
}



int samRead::printAln(FILE *nf, std::string &chrSeq){
  map<int,int> mc;
  this->getMapOfCoords(mc);
  return ( this->printAln(nf, mc, chrSeq) );
}


int samRead::printAln(FILE *nf, std::map<int,int> &m, std::string &chrSeq){

  char bRead, bChr;

  map<int,int>::iterator it;
  for(it=m.begin() ; it!=m.end() ; it++){
    int posInRead = it->first;
    int posOnChr = it->second;
    if(posInRead != -1){ bRead = this->seq_[posInRead]; } else { bRead = '-'; }
    if( posOnChr != -1){ bChr = chrSeq[posOnChr]; } else { bChr = '-'; }
    fprintf(nf, "[%d - %d]\t%c:%c\n", posInRead, posOnChr, bRead, bChr);
  }
  return 0;
}


int samRead::tryToFillDeletionLeft(map<int, int> &newM, int indelPos, int indelSize, GFFtranscript &gt, string &chrSeq){

  int iPos;
  int posInRead = indelPos;
  int readSize = newM.size();
  int chrLg = chrSeq.size();
  int transcriptLg = gt.getSplicedSize();

  int shift = 1;
  if(gt.isOnPlusStrand()==false){
    shift = -1;
  }

  int posGenomic = newM[indelPos];
  int lastPosGenomic = posGenomic;

  int posInTranscript = gt.genomicToSplicedTranscript(posGenomic) + shift;
  posGenomic = gt.splicedToGenomic(posInTranscript);
  posInRead++;

  if(posGenomic<0){ // If I'm outside an annotated gene
    posGenomic = lastPosGenomic + 1;
    lastPosGenomic = posGenomic;
  }

  char baseGenomic = 'N';
  if( posGenomic>=0 && posGenomic<chrLg ){
    baseGenomic = chrSeq[posGenomic];
  }

  int nbMoved = 0;
  while( posInRead<readSize && baseGenomic == this->seq_[posInRead] ){

    nbMoved++;
    newM[posInRead] = posGenomic;

    posInRead++;
    posInTranscript += shift;
    posGenomic = gt.splicedToGenomic(posInTranscript);

    if(posGenomic<0){ // If I'm outside an annotated gene
      posGenomic = lastPosGenomic + 1;
      lastPosGenomic = posGenomic;
    }

    if( posGenomic>=0 && posGenomic<chrLg ){
      baseGenomic = chrSeq[posGenomic];
    } else {
      baseGenomic = 'N';
    }

    
  }

  return nbMoved;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//                                    UTILITY FUNCTIONS                                         //
//////////////////////////////////////////////////////////////////////////////////////////////////


void samRead_vectoriseCigar(vector< pair<char, int> > &v, const string &s){
  char ctmp[8];

  int cigarPos, subPos;

  int lg = s.size();

  for(cigarPos=0, subPos=0 ; cigarPos<lg ; cigarPos++){
    char c = s[cigarPos];
    if(c<'0' || c>'9'){
      ctmp[subPos] = '\0';
      int nb = atoi(ctmp);
      pair<char, int> pp(c, nb);
      v.push_back(pp);
      subPos = 0;
    } else {
      ctmp[subPos] = s[cigarPos];
      subPos++;
    }
  }
}



void myLocal_samRead_printMap(FILE *nf, map<int,int> & m){
  map<int,int>::iterator it;
  for(it=m.begin() ; it!=m.end() ; it++){
    fprintf(nf, "%d -> %d\n", it->first, it->second);
  }
}


bool samRead::isOnIntron(void){
  string s = this->getCigar();
  if(s.find('N')!=string::npos)
    return true;
  return false;
}

#endif
