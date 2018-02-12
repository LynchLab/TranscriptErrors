/*
CS_consensus.cxx
 */


#include "naiveRepeatSearch.hxx"
#include "CS_consensus.hxx"


using namespace std;

//void myLocal_printMap(FILE *nf, map<int,int> &m, string &cs, string &chrSeq);

CS_consensus::~CS_consensus(void){
  this->clean();
}


void CS_consensus::clean(void){
  this->mc_.clear();
}


/**
   Default constructor, does not do anything.
*/
CS_consensus::CS_consensus(void){
    this->consensus_ = "";
}


/**
   Constructs a consensus sequence from a read
   @param read the read sequence
   @param qual quality of the read
   @param minRepeatSize minimum size for repeat to search
   @param maxMM maximum number of mismatches between the priming region (position 1 to minRepeatSize in read) and the first repeat occurence.
 */
CS_consensus::CS_consensus(const std::string &read, const std::string &qual, const int minRepeatSize, const int maxMM=2, const double minID=0.9){
  this->make(read, qual, minRepeatSize, maxMM, minID);
}


/**
   Constructs a consensus sequence from a pair of reads
   @param leftRead the left-read sequence
   @param leftQual quality of the left read
   @param rightRead the left-read sequence
   @param rightQual quality of the right read
   @param minRepeatSize minimum size for repeat to search
   @param maxMM maximum number of mismatches between the priming region (position 1 to minRepeatSize in read) and the first repeat occurence.
 */
CS_consensus::CS_consensus(const std::string &leftRead, const std::string &leftQual, const std::string &rightRead, const std::string &rightQual, const int minRepeatSize, const int maxMM=2, const double minID=0.9){
  this->make(leftRead, leftQual, rightRead, rightQual, minRepeatSize, maxMM, minID);
}


/**
   Constructs a consensus from a read (with the consensus sequence already known)
   @param read the read sequence
   @param qual quality of the read
   @param csSeq the consensus sequence
   @param minID the minimum fraction of identical nucleotides for a subsequence of the read to be considered as a unit of repeat
 */
CS_consensus::CS_consensus(const std::string &read, const std::string &qual, const std::string csSeq, const double minID=0.9){
  this->make(read, qual, csSeq, minID);
}


/**
   Constructs a consensus from a pair of reads (with the consensus sequence already known)
   @param leftRead the left-read sequence
   @param leftQual quality of the left-read
   @param rightRead the right-read sequence
   @param rightQual quality of the right-read
   @param csSeq the consensus sequence
   @param minID the minimum fraction of identical nucleotides for a subsequence of the read to be considered as a unit of repeat
 */
CS_consensus::CS_consensus(const std::string &leftRead, const std::string &leftQual, const std::string &rightRead, const std::string &rightQual, const std::string csSeq, const double minID=0.9){
  this->make(leftRead, leftQual, rightRead, rightQual, csSeq, minID);
}




/**
   The function that is called by the constructor for a single read.
   Requires to find at least one full repeat (=2 occurrences) with id>=minID
   @param Read the read sequence
   @param Qual the read quality scores
   @param minRepeatSize minimum size for repeat to search
   @param maxMM maximum number of mismatches between the priming region (position 1 to minRepeatSize in read) and the first repeat occurence.
   @param minID the minimum fraction of identical nucleotides for a subsequence of the read to be considered as a unit of repeat
*/
bool CS_consensus::make(const std::string &Read, const std::string &Qual, const int minRepeatSize, const int maxMM, const double minID){

  //int pos = findFirstRepeatNaive(Read, minRepeatSize, maxMM);

  int nbMM = -1;
  // First, I search a repeat of the first 'minRepeatSize' bases
  string stmp = Read.substr(0, minRepeatSize);
  int pos = findFirstOccurrenceNaive(stmp, Read, maxMM, &nbMM, minRepeatSize);

  if( pos > 0 ){
    int repeatSize = pos;
    string firstOccurrence = Read.substr(0, repeatSize);

    int maxMMtot = repeatSize * (1.0-minID);

    // Searching for the second occurrence
    if( findFirstOccurrenceNaive(firstOccurrence, Read, maxMMtot, &nbMM, repeatSize) > 0 ){
      // Valid repeat
      this->posStartFirst_.first = pos;
      this->posStartFirst_.second = -1;
      this->updateMap(Read, Qual, 0, repeatSize);
      this->updateConsensus();
      return true;
    }

  }

  this->mc_.clear();
  this->consensus_ = "";
  this->posStartFirst_.first = -1;
  this->posStartFirst_.second = -1;
  return false;

}



/**
   The function that is called by the constructor for a pair of reads.
   This function starts by searching a consensus in the left read. After that, it searches for an identical consensus in the right read and builds the overall consensus only if consensus is identical between left and right reads.
   @param leftRead the left-read sequence
   @param leftQual the left-read quality scores
   @param rightRead the right-read sequence
   @param rightQual the right-read quality scores
   @param minRepeatSize minimum size for repeat to search
   @param maxMM maximum number of mismatches between the priming region (position 1 to minRepeatSize in read) and the first repeat occurence.
   @param minID the minimum fraction of identical nucleotides for a subsequence of the read to be considered as a unit of repeat
 */
bool CS_consensus::make(const std::string &leftRead, const std::string &leftQual, const std::string &rightRead, const std::string &rightQual, const int minRepeatSize, const int maxMM, const double minID){

  // First, I try to build a consensus from the left read.
  if( this->make(leftRead, leftQual, minRepeatSize, maxMM, minID) == false )
    return false;

  //fprintf(stderr, "left-read consensus:\t%s\n", this->getConsensus().c_str());

  // Searching for the consensus (found in the left read) in the right read
  // Requirement = at least one occurrence
  string leftCs = this->getConsensus();
  int nbMM = -1;
  int maxMMtot = (1.0 - minID) * leftCs.size();
  int pos = findFirstOccurrenceNaive(leftCs, rightRead, maxMMtot, &nbMM, 0);

  //  fprintf(stderr, "DEBOG:\n -Find: '%s'\n -In:'%s'\n -pos:%d (maxMM:%d)\n", leftCs.c_str(), rightRead.c_str(), pos, maxMM);

  if( pos > 0 ){
    int repeatSize = leftCs.size();

    if( pos >= repeatSize ){
      pos = pos % repeatSize;
    }

    this->posStartFirst_.second = pos;
    this->updateMap(rightRead, rightQual, pos, repeatSize);
    this->updateConsensus();
    return true;
  } else {
    this->posStartFirst_.first = -1;
    this->posStartFirst_.second = -1;
    this->mc_.clear();
    this->consensus_ = "";
    return false;
  }

  this->posStartFirst_.first = -1;
  this->posStartFirst_.second = -1;

  return false;
}



bool CS_consensus::make(const string &read, const string &qual, const string &csSeq, const double minID=0.9){

  int nbMM = -1;
  int maxMM = (1.0 - minID) * csSeq.size();
  int pos = findFirstOccurrenceNaive(csSeq, read, maxMM, &nbMM, 0);

  if( pos<0 ){ return false; }
  if( pos>= csSeq.size() ){
    pos = pos % csSeq.size();
  }
  this->updateMap(read, qual, pos, csSeq.size());
  this->updateConsensus();
  return true;

}


bool CS_consensus::make(const std::string &leftRead, const std::string &leftQual, const std::string &rightRead, const std::string &rightQual, const std::string &csSeq, const double minID=0.9){

  if( this->make(leftRead, leftQual, csSeq, minID) == false ){
    this->mc_.clear();
    this->consensus_ = "";
    return false;
  }

  if( this->make(rightRead, rightQual, csSeq, minID) == false ){
    this->mc_.clear();
    this->consensus_ = "";
    return false;
  }

  return true;

}



bool CS_consensus::make(const std::string &readSeq, const std::string &readQual, const std::string &mateSeq, const std::string &mateQual, const int csLg, const int breakPos, const int csPosInMate){
  this->updateMap(readSeq, readQual, breakPos, csLg);
  this->updateConsensus();
  this->updateMap(mateSeq, mateQual, (csPosInMate+breakPos)%csLg, csLg);
  this->updateConsensus();
}


bool CS_consensus::make(const std::string &readSeq, const std::string &readQual, const int csLg, const int breakPos){
  this->updateMap(readSeq, readQual, breakPos, csLg);
  this->updateConsensus();
}



/**
   Call this function after finding the first occurence of the repeat to fill in the map
   @parama seq The read sequence
   @param qual The read quality
   @param posFirstRepeat The position of the first repeat in the read
   @param repeatSize The size of the repeat
 */
void CS_consensus::updateMap(const std::string &seq, const std::string &qual, const int posFirstRepeat, const int repeatSize){

  int posInConsensus;

  //fprintf(stderr, "updateMap debog:\nseq: '%s'\nposFirstRepeat:%d\nrepeatSize:%d\n", seq.c_str(), posFirstRepeat, repeatSize);

  int seqSize = seq.size();
  for(int posInSeq=0 ; posInSeq<seqSize ; posInSeq++){
    posInConsensus = (posInSeq - posFirstRepeat) % repeatSize;
    posInConsensus = ((repeatSize - posFirstRepeat) + posInSeq) % repeatSize;
    this->mc_[posInConsensus].first.push_back(seq[posInSeq]);
    this->mc_[posInConsensus].second.push_back(qual[posInSeq]);
  }

}




/**
   Updates the consensus sequence based on the list of base calls and quality scores at each position.
 */
void CS_consensus::updateConsensus(void){

  int id, ib, depth, posInConsensus;

  this->consensus_ = "";

  map<int, pair<string,string> >::iterator it;
  for(it=this->mc_.begin() ; it!=this->mc_.end() ; it++){
    posInConsensus = it->first;

    // baseAtPos and qualAtPos contain all the calls and qualities at a given position in a consensus, concatenated in strings
    string basesAtPos = (it->second).first;
    string qualAtPos = (it->second).second;

    char tabBases[5] = {'A', 'T', 'C', 'G', 'N'};
    int tabQual[5] = {0, 0, 0, 0, 0};

    depth = basesAtPos.size();
    for(id=0 ; id<depth ; id++){
      for(ib=0 ; ib<4 && tabBases[ib]!=basesAtPos[id]; ib++)
	;
      tabQual[ib] += ( (qualAtPos[id] - '"') + 1); // The illumina lowest score (1) corresponds to ASCII character double quote ("). Here I transform the ASCII character into a quality score (from 1 to 41)
    }

    // Searching for base with the higest cumulative quality score
    int iMaxQual = 4;
    for(ib=0 ; ib<4 ; ib++){
      if(tabQual[ib]>tabQual[iMaxQual]){
	iMaxQual = ib;
      }
    }

    this->consensus_.push_back(tabBases[iMaxQual]);

  }

}


/**
   Computes the quality scores of the consensus sequence.
 */
string CS_consensus::getCsQuality(const int minQual, const int maxQual){

  int lg = this->consensus_.size();
  string csQual(lg, minQual);

  map<int, pair<string,string> >::iterator itm;
  for(itm=this->mc_.begin() ; itm!=this->mc_.end() ; itm++){
    int pos = itm->first;
    int scoreCs = 0;
    int scoreNonCs = 0;

    char cs = this->consensus_[pos];
    string sObs = itm->second.first;
    string sQual = itm->second.second;
    int depth = sObs.size();

    for(int iDepth=0 ; iDepth<depth ; iDepth++){
      if( sObs[iDepth]==cs ){
	scoreCs += sQual[iDepth];
      } else {
	scoreNonCs += sQual[iDepth];
      }
    }

    int scoreLocal = scoreCs - scoreNonCs;
    if( scoreLocal < minQual ){ scoreLocal = minQual; }
    if( scoreLocal > maxQual ){ scoreLocal = maxQual; }
    csQual[pos] = scoreLocal;

  }

  return csQual;
}




/**
   Computes the confidence score associated to each base in the consensus sequence. The confidence scores are stored in the vector of integer (passed as reference).
   Each integer represents the Phred score, which can be converted to a probability with the forumla: P = 10^-Q/10 (Q = integer)
   
 */
int CS_consensus::getCsProb(vector<int> &vQual){

  int i;

  int lg = this->consensus_.size();

  vQual.clear();
  for(i=0 ; i<lg ; i++){
    vQual.push_back(0);
  }

  map<int, pair<string,string> >::iterator itm;
  for(itm=this->mc_.begin() ; itm!=this->mc_.end() ; itm++){
    int pos = itm->first;
    int scoreCs = 0;
    int scoreNonCs = 0;

    char cs = this->consensus_[pos];
    string sObs = itm->second.first;
    string sQual = itm->second.second;
    int depth = sObs.size();

    for(int iDepth=0 ; iDepth<depth ; iDepth++){
      if( sObs[iDepth]==cs ){
	scoreCs += sQual[iDepth];
      } else {
	scoreNonCs += sQual[iDepth];
      }
    }

    int scoreLocal = scoreCs - scoreNonCs;
    //if( scoreLocal < minQual ){ scoreLocal = minQual; }
    //if( scoreLocal > maxQual ){ scoreLocal = maxQual; }
    //csQual[pos] = scoreLocal;

    if(pos<0 || pos>=lg){
      fprintf(stderr, "PROBLEM: pos=%d - lg=%d\n", pos, lg);
      fprintf(stderr, "'%s'\n", this->consensus_.c_str());
      fflush(stderr);
      return -1;
    } else {
      vQual[pos] = scoreLocal;
    }

  }
  return 0;
}


/**
   void CS_consensus::getCsProb(vector<int> &vQual, vector<int> &vDepth, vector<int> vDivScores){
   Computes the confidence score associated to each base in the consensus sequence. The confidence scores are stored in the vector of integer (passed as reference).
   Each integer represents the Phred score, which can be converted to a probability with the forumla: P = 10^-Q/10 (Q = integer)
   Also fills in two additional vectors:
    - vDepth: to record the depth of the consensus (= how many repeats) at each position.
    - vDivScores: the sum of divergent base calls scores at each position.
   
 */
int CS_consensus::getCsProb(vector<int> &vQual, vector<int> &vDepth, vector<int> &vDivScores){

  int i;

  int lg = this->consensus_.size();

  vQual.clear();
  vDepth.clear();
  vDivScores.clear();
  for(i=0 ; i<lg ; i++){
    vQual.push_back(0);
    vDepth.push_back(0);
    vDivScores.push_back(0);
  }

  map<int, pair<string,string> >::iterator itm;
  for(itm=this->mc_.begin() ; itm!=this->mc_.end() ; itm++){
    int pos = itm->first;
    int scoreCs = 0;
    int scoreNonCs = 0;

    char cs = this->consensus_[pos];
    string sObs = itm->second.first;
    string sQual = itm->second.second;
    int depth = sObs.size();

    for(int iDepth=0 ; iDepth<depth ; iDepth++){
      if( sObs[iDepth]==cs ){
	scoreCs += sQual[iDepth];
      } else {
	scoreNonCs += sQual[iDepth];
      }
    }

    vDepth[pos] = depth;
    vDivScores[pos] = scoreNonCs;

    int scoreLocal = scoreCs - scoreNonCs;
    //if( scoreLocal < minQual ){ scoreLocal = minQual; }
    //if( scoreLocal > maxQual ){ scoreLocal = maxQual; }
    //csQual[pos] = scoreLocal;

    if(pos<0 || pos>=lg){
      fprintf(stderr, "PROBLEM: pos=%d - lg=%d\n", pos, lg);
      fprintf(stderr, "'%s'\n", this->consensus_.c_str());
      fflush(stderr);
      return -1;
    } else {
      vQual[pos] = scoreLocal;
    }

  }

  return 0;
}








/**
   Prints the consensus in a vertical fashion.
   @param nf Pointer to the output stream (opened file, typically stdout or stderr)
   @param printQual True -> print quality scores / False -> DO not print quality scores
 */
void CS_consensus::prettyPrint(FILE *nf, const bool printQual=false){

  int i, j, localDepth;
  int lg = this->consensus_.size();
  int maxDepth = -1;

  for(i=0 ; i<lg ; i++){
    localDepth = this->mc_[i].first.size();
    if( localDepth > maxDepth ){
      maxDepth = localDepth;
    }
  }

  //fprintf(nf, "%s\n", this->consensus_.c_str());

  for(j=0 ; j<maxDepth ; j++){
    for(i=0 ; i<lg ; i++){
      localDepth = this->mc_[i].first.size();
      if( j>=localDepth ){
	fprintf(nf, "%c", ' ');
      } else {
	fprintf(nf, "%c", this->mc_[i].first[j]);
      }
    }
    fprintf(nf, "\n");
  }

  if( printQual==true ){

    for(j=0 ; j<maxDepth ; j++){
      for(i=0 ; i<lg ; i++){
	localDepth = this->mc_[i].second.size();
	if( j>=localDepth ){
	  fprintf(nf, "%c", ' ');
	} else {
	  fprintf(nf, "%c", this->mc_[i].second[j]);
	}
      }
      fprintf(nf, "\n");
    }

  }

}


int CS_consensus::getPosStartInRead(void){
  return this->posStartFirst_.first;
}


int CS_consensus::getPosStartInMate(void){
  return this->posStartFirst_.second;
}


void CS_consensus::reverseCpt(const char *tab_cpt){

  map<int, pair<string, string> > newMap;

  map<int, pair<string, string> >::reverse_iterator itm;
  int pos = 0;

  for(itm=this->mc_.rbegin(), pos=0 ; itm!=this->mc_.rend() ; itm++, pos++){
    string sCalls = (itm->second).first;
    string sQuals = (itm->second).second;

    int depth = sCalls.size();
    for(int iDepth=0 ; iDepth<depth ; iDepth++){
      sCalls[iDepth] = tab_cpt[sCalls[iDepth]];
    }

    pair<string,string> pp(sCalls, sQuals);
    newMap[pos] = pp;
  }

  this->mc_ = newMap;

  this->updateConsensus();

}





int CS_consensus::refine(int *offset, samRead &sr, vector<GFFtranscript> &vt, map<string,string> &mg){

  int i, j;
  //int minEdit = sr.getEditDistance();
  int initialIndel = sr.getNbIndels();
  int initialMM = sr.getNbMM();
  int initialEditDist = initialIndel + initialMM;

  int minIndel = -1;
  int minMM = -1;
  int minEdit = -1;
  int offsetFromMinEdit = -1;


  map<int,int> mc;
  if( sr.getMapOfCoords(mc) != 0 ){
    fprintf(stderr, "WARNING: FAILLED TO OBTAIN map of coordinates for '%s'\n", sr.getQueryName().c_str());
    return -1;
  }

  //myLocal_samRead_printMap(stderr, mc);

  int nbt = vt.size();
  if(nbt<1){
    fprintf(stderr, "%s does not map on an annotated gene.\n", sr.getQueryName().c_str());
    return -1;
  }
  for(i=0 ; i<nbt ; i++){
    GFFtranscript gt = vt[i];
    if(gt.chr_ != sr.getRefName()){
      fprintf(stderr, 
	      "WARNING: read chromosome != GFFtranscript chromosome. SHould never happen. readID:'%s' - transcriptID:'%s'\n",
	      sr.getQueryName().c_str(),
	      gt.id_.c_str()
	      );
    }
    int tmpOffset;
    int editDist = sr.refine(mc, mg[gt.chr_], gt, &tmpOffset);

    //fprintf(stderr, "%s\n", gt.id_.c_str());

    if(editDist>-1){
      
      if(editDist>initialEditDist){
	//fprintf(stderr, "NO REFINE: %s\n", sr.getQueryName().c_str());
	;
	// ADD a continue statement?
      }
      
      // The third condition in the next if
      if( 
	 (editDist<initialEditDist && minEdit==-1) || // First try
	 editDist<minEdit || // Improvement in the editDist
	 (editDist==initialEditDist && (editDist<minEdit || minEdit==-1) && initialIndel>0 )  //Total editDistance does not change but move from indel to mismatch
	  ){
	minEdit = editDist;
	offsetFromMinEdit = tmpOffset;
      }
    }
    
  }


  if( minEdit>=0 &&( minEdit<initialEditDist || (minEdit==initialEditDist && initialIndel>0) ) ){

    map<int, pair<string,string> > mcs;

    this->updateConsensus();
    int csLg = this->consensus_.size();
    
    // Here, update the consensus map
    if( offsetFromMinEdit==0 ){
      fprintf(stderr, "ERROR: OFFSET==0. ReadID:'%s'\n", sr.getQueryName().c_str());
      return -1;
    }
    
    if( offsetFromMinEdit>0 ){
      for(i=offsetFromMinEdit ; i<csLg ; i++){
	mcs[i-offsetFromMinEdit] = this->mc_[i];
      }
      for(j=0 ; j<offsetFromMinEdit ; j++, i++){
	mcs[i-offsetFromMinEdit] = this->mc_[j];
      }
    } else {
      int posBreak = csLg - (-offsetFromMinEdit);
      for(i=0, j=posBreak ; i<(-offsetFromMinEdit) ; i++, j++){
	mcs[i] = this->mc_[j];
      }
      for(j=0 ; j<posBreak ; j++, i++){
	mcs[i] = this->mc_[j];
      }
    }
    
    this->mc_.clear();
    this->mc_ = mcs;
    this->updateConsensus();
    
    //fprintf(stderr, "New consensus: %s\n", this->consensus_.c_str());
    
    (*offset) = offsetFromMinEdit;
    return minEdit;

  }

  return -2;
  
}


/**
   Computes the identity percentage between the consensus sequence and each of the repeats.

 */
double CS_consensus::getPrctId(void){

  int i;
  int nbTot = 0;
  int nbIdent = 0;
  int nbDiff = 0;

  int csLg = this->consensus_.size();

  map<int, pair<string,string> >::iterator itm;
  for(itm=this->mc_.begin() ; itm!=this->mc_.end() ; itm++){
    int pos = itm->first;
    string sTmp = (itm->second).first;
    int localDepth = sTmp.size();
    for(i=0 ; i<localDepth ; i++){
      if(sTmp[i]==this->consensus_[pos]){
	nbIdent++;
      } else {
	nbDiff++;
      }
      nbTot++;
    }
  }

  double id = (double)((double)nbIdent/(double)nbTot);

  return id;
}

/**
   Computes the how far homopolymer tracks run from each end of the consensus sequence.
 */
void CS_consensus::setLimitHomo(void){
  int i;
  int csLg = this->consensus_.size();
  for(i=1 ; i<csLg && this->consensus_[i] == this->consensus_[0] ; i++)
    ;
  this->limitHomo_.first = i;

  int lastPos = csLg - 1;
  for(i=lastPos-1 ; i>=0 && this->consensus_[i] == this->consensus_[lastPos] ; i--)
    ;
  this->limitHomo_.second = csLg - i;
}


void CS_consensus::setLimitHomo(int left, int right){
  this->limitHomo_.first = 0;
  this->limitHomo_.second = 0;
}


/**
   Gets the size of homopolymers that are to be exclude at both ends of the consensus sequence.
 */
void CS_consensus::getHomoSizes(int *left, int *right){
  int csLg = this->consensus_.size();
  (*left) = this->limitHomo_.first;
  (*right) = this->limitHomo_.second;
}

void myLocal_printMap(FILE *nf, map<int,int> &m, string &cs, string &chrSeq){

  map<int, int>::iterator it;
  for(it=m.begin() ; it!=m.end() ; it++){
    int posInCs = it->first;
    int posGenomic = it->second;
    char csBase = cs[posInCs];
    char genomicBase = chrSeq[posGenomic];
    fprintf(nf, "%-2d\t%c %c\n", posInCs, csBase, genomicBase);
  }

}

/**
   Returns the length of the consensus sequence.
 */
int CS_consensus::getSize(void){
  return this->consensus_.size();
}
