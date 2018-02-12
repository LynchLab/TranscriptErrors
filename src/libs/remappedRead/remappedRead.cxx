/*

remappedRead.cxx

 */


#include <iostream>
#include <fstream>
#include <sstream>


#include "remappedRead.hxx"



using namespace std;


remappedRead::remappedRead(void){
  ;
}


remappedRead::remappedRead(string s, bool includeSequence, bool computeSeqWithIndels, bool includeQuality){

  string stmp, spair, sstart, slg;
  long start, lg;

  if( includeSequence==false && (includeQuality==true || computeSeqWithIndels==true) ){
    if(includeQuality==true){
      fprintf(stderr, "ERROR: quality cannot be included into the remappedRead without the sequence...\n");
      exit(1);
    }
    if(computeSeqWithIndels==true){
      fprintf(stderr, "ERROR: impossible to compute sequence with indels without the sequence...\n");
      exit(1);
    }
  }

  istringstream iss(s);

  getline(iss, this->qName_, '\t');
  getline(iss, this->strand_, '\t');
  getline(iss, this->rName_, '\t');

  getline(iss, this->transcriptID_, '\t');

  // List of (start:lg,) pairs of mapping coordinates on the genome
  getline(iss, stmp, '\t');
  istringstream issMapCoords(stmp);
  while( getline(issMapCoords, spair, ',') ){
    istringstream iss2(spair);
    getline(iss2, sstart, ':');
    getline(iss2, slg, ':');
    pair<long, long> pp;
    pp.first = atol(sstart.c_str());
    pp.second = atol(slg.c_str());
    this->vmaps_.push_back(pp);
  }


  getline(iss, stmp, '\t');
  if(stmp.find("deletion:")!=string::npos){ stmp.erase(0, 9); }
  istringstream issDeletions(stmp);
  while( getline(issDeletions, spair, ',') ){
    istringstream iss2(spair);
    getline(iss2, sstart, ':');
    getline(iss2, slg, ':');
    pair<long, long> pp;
    pp.first = atol(sstart.c_str());
    pp.second = atol(slg.c_str());
    this->vDeletions_.push_back(pp);
  }


  getline(iss, stmp, '\t');
  if(stmp.find("insertion:")!=string::npos){ stmp.erase(0, 10); }
  istringstream issInsertions(stmp);
  while( getline(issInsertions, spair, ',') ){
    istringstream iss2(spair);
    getline(iss2, sstart, ':');
    getline(iss2, slg, ':');
    pair<long, long> pp;
    pp.first = atol(sstart.c_str());
    pp.second = atol(slg.c_str());
    this->vInsertions_.push_back(pp);
  }

  if( includeSequence == true ){
    getline(iss, this->seq_, '\t');

    if( computeSeqWithIndels == true ){
      this->seqWithIndels_ = this->getSeqWithIndels();
    }

    if( includeQuality == true ){
      getline(iss, this->squal_, '\t');
    }
  }

}


string remappedRead::getSeqWithIndels(bool qual){

  string seqWithIndels = this->seq_;
  if(qual==true){
    seqWithIndels = this->squal_;
  }

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Now I will modify the sequence of the read so as to introduce the indels. WARNING: this code works only if there is only one insertion/deletion per read !!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  for(int i=0 ; i<vDeletions_.size() ; i++){

    if( i>0 ){
      fprintf(stderr, "ERROR: multiple deletions.\n");
      return "";
    }

    pair<long, long> pp = vDeletions_[i];
    long posDeletion = pp.first;
    long sizeDeletion = pp.second;
    long posDeletionInRead = this->genomeToReadCoords(posDeletion);
    string si(sizeDeletion, '-');

    if( posDeletionInRead<0 ){
      fprintf(stderr, 
	      "Error while calling this->genomeToReadCoords(%d) for %s\n", 
	      posDeletion,
	      this->qName_.c_str()
	      );
      return "";
    }

    seqWithIndels.insert(posDeletionInRead, si);
  }
  

  for(int i=0 ; i<vInsertions_.size() ; i++){

    if( i>0 ){
      fprintf(stderr, "ERROR: multiple insertions.\n");
      return "";
    }

    pair<long, long> pp = vInsertions_[i];
    long posInsertion = pp.first;
    long sizeInsertion = pp.second;
    long posInsertionInRead = this->genomeToReadCoords(posInsertion);

    if( posInsertionInRead<0 ){
      fprintf(stderr, 
	      "Error while calling this->genomeToReadCoords(%d) for %s\n", 
	      posInsertion,
	      this->qName_.c_str()
	      );
	      
      return "";
    }

    seqWithIndels.erase(posInsertionInRead, sizeInsertion);
  }

  return seqWithIndels;

}




long remappedRead::genomeToReadCoords(long pos){

  long posInRead;

  if( this->coversPosition(pos) == false ){
    return -1;
  }

  if(this->strand_=="+"){
    long nbToSubstract = this->vmaps_[0].first;
    for(int i=1 ; i < this->vmaps_.size() && pos >= this->vmaps_[i].first ; i++){
      long distance = vmaps_[i].first - (vmaps_[i-1].first + vmaps_[i-1].second);
      nbToSubstract += distance;
    }

    posInRead = pos - nbToSubstract;

  } else {
    long nbToSubstractFrom = this->vmaps_[0].first;
    for(int i=1 ; i < this->vmaps_.size() && pos <= this->vmaps_[i].first ; i++){
      long distance = (this->vmaps_[i-1].first - this->vmaps_[i-1].second) - this->vmaps_[i].first;
      nbToSubstractFrom -= distance;
    }

    posInRead = nbToSubstractFrom - pos;
  }

  return posInRead;
}


void remappedRead::quickDisplay(FILE *nf){

  fprintf(nf, "%s\t%s\t%d\n", 
	  this->qName_.c_str(),
	  this->rName_.c_str(),
	  this->vmaps_[0].first
	  );

}


bool remappedRead::coversPosition(long pos){

  int i;
  int nbBlocs = this->vmaps_.size();
  bool cpt = false;
  if(this->strand_=="-"){ cpt = true; }


  for(i=0 ; i<nbBlocs ; i++){
    pair<long, long> pp = this->vmaps_[i];
    long start, end;
    if(cpt==true){
      start = pp.first - pp.second;
      end = pp.first;
    } else {
      start = pp.first;
      end = start + pp.second;
    }

    if(start<pos && end>pos){
      return true;
    }

  }

  return false;
}



void remappedRead::fullPrint(FILE *nf, bool printQuality){

  int i;

  fprintf(nf, "%s\t%s\t%s\t%s\t", 
	  this->qName_.c_str(), 
	  this->strand_.c_str(), 
	  this->rName_.c_str(),
	  this->transcriptID_.c_str()
	  );

  for(i=0 ; i<this->vmaps_.size() ; i++){
    pair<long, long> pp = this->vmaps_[i];
    fprintf(nf, "%d:%d,", pp.first, pp.second);
  }

  // Printing out the list of deletions:
  fprintf(nf, "\t");
  for(i=0 ; i<this->vDeletions_.size() ; i++){
    pair<long, long> pp = this->vDeletions_[i];
    fprintf(nf, "deletion:%d:%d,", pp.first, pp.second);
  }

  // Printing out the list of insertions:
  fprintf(nf, "\t");
  for(i=0 ; i<this->vInsertions_.size() ; i++){
    pair<long, long> pp = this->vInsertions_[i];
    fprintf(nf, "insertion:%d:%d,", pp.first, pp.second);
  }

  fprintf(nf, "\t%s", this->seq_.c_str());

  if(printQuality==true){
    fprintf(nf, "\t%s", this->squal_.c_str());
  }

}



bool remappedRead::asInsertion(pair<long, long> & pi){

  if( find(this->vInsertions_.begin() , this->vInsertions_.end() , pi) == this->vInsertions_.end() ){
    return false;
  }

  return true;
}


bool remappedRead::asDeletion(pair<long, long> & pi){

  if( find(this->vDeletions_.begin() , this->vDeletions_.end() , pi) == this->vDeletions_.end() ){
    return false;
  }

  return true;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
// UTILITY FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////


