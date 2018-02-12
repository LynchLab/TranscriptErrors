#include "naiveRepeatSearch.hxx"
#include "UTILS_FASTA/utils_fasta.hxx"


using namespace std;


int findFirstRepeatNaive(const string &s, const int minRepeatSize, const int maxMM){

  int mpos, spos, nbMM;

  string motif = s.substr(0, minRepeatSize);
  int maxStart = (s.size() - minRepeatSize) + 1; // The right-most position to start searching for a repeat


  for(spos=minRepeatSize ; spos<=maxStart ; spos++){
    for(mpos=0, nbMM=0 ; nbMM<=maxMM && mpos<minRepeatSize ; mpos++){
      if( motif[mpos]!=s[spos+mpos] ){
	nbMM++;
      }
    }
    
    if(nbMM<=maxMM && mpos==minRepeatSize){
      string firstRepeat = s.substr(minRepeatSize, (spos+mpos)-minRepeatSize);
      int nbMMSecondRepeat = -1;
      int posSecondRepeat = findFirstOccurrenceNaive(firstRepeat, s, maxMM, &nbMMSecondRepeat, (spos+mpos));
      if( posSecondRepeat > 0 ){ return spos; }
    }

  }
  return -1;
}


int findSeqWithMM(const string &seqToFind, const string &seq, const int maxMM){

  int maxStart = seq.size() - seqToFind.size();
  int lgf = seqToFind.size();
  int i, j, nbMM;

  for(i=0 ; i<maxStart ; i++){
    for(j=0, nbMM=0 ; j<lgf && nbMM<=maxMM ; j++){
      if( seq[i+j]!=seqToFind[j] ){
	nbMM++;
      }
    }

    if(nbMM<=maxMM && j==lgf)
      return i;

  }

  return -1;
}


/** 
    Finds the first occurence of a string within a string, allowing for mismatches and different start position for search.
    @param sToFind The string to search for
    @param sToSearch The string into which the search is performed
    @param maxMM Maximum number of mismatches between the string to find and the hit
    @nbMM Pointer to the variable that will be filled with the number of mismatches
    @startSearchFromPos Start searching from this position
    @return The position of the first occurence found.
*/
int findFirstOccurrenceNaive(const std::string &sToFind, const std::string &sToSearch, const int maxMM, int *nbMM, const int startSearchFromPos){

  int startFrom = startSearchFromPos;
  if(startSearchFromPos<0){
    startFrom = 0;
  }
  
  int searchSize = sToSearch.size();
  int findSize = sToFind.size();

  int posInSearch, posInFind, nbMM_local;

  int maxPosInSearch = (searchSize - findSize) + 1;

  for(posInSearch = startFrom ; posInSearch <= maxPosInSearch ; posInSearch++){
    nbMM_local = 0;
    posInFind = 0;
    while( posInFind<findSize && nbMM_local<=maxMM ){
      char ccSearch = sToSearch[posInSearch + posInFind];
      char ccFind = sToFind[posInFind];
      if( ccSearch != ccFind ){
	nbMM_local++;
      }
      posInFind++;
    }

    if( posInFind==findSize && nbMM_local<=maxMM ){
      (*nbMM) = nbMM_local;
      return posInSearch;
    }

  }

  return -1;
}






/**
   Searches for a previously found consensus sequence in a new read (typically used to search for left-read-derived consensus in the right-read sequence)
   @param consensus The previously derived consensus sequence
   @param read the read sequence
   @param qual quality of the read
   @param maxMM maximum number of mismatches between the previously defined consensus sequence and the best match in the new sequence
 */
bool findConsensus(const std::string &consensus, const std::string &read, const std::string qual, const int maxMM){

  int csSize = consensus.size();
  int seqSize = read.size();
  int nbMM;

  int bestConsensusPos = -1;
  int bestMM = maxMM + 1;
  int lastConsensusFound = findFirstOccurrenceNaive(consensus, read, maxMM, &nbMM);

  while( lastConsensusFound >= 0 ){
    if(nbMM<bestMM){
      bestMM = nbMM;
      bestConsensusPos = lastConsensusFound;
    }
    lastConsensusFound = findFirstOccurrenceNaive(consensus, read, maxMM, &nbMM, lastConsensusFound+1);
  }

  if( bestConsensusPos >= 0 ){
    return bestConsensusPos;
  }

  return false;

}

