/*

GFFtranscript.cxx

 */

#include <algorithm>

#include "UTILS_FASTA/utils_fasta.hxx"
#include "GFFtranscripts.hxx"

using namespace std;

GFFtranscript::GFFtranscript(void){

  this->posStart_ = -1;
  this->posStop_ = -1;
  
}



GFFtranscript::GFFtranscript(vector<GFFentry> & vEntries, string sTAG, int shift, long posStart, long posStop){
  map<string, string> mg;
  this->init(vEntries, sTAG, shift, mg, NULL, posStart, posStop);
}


// vEntries: vector of GFFentry all from the same sTAG (and ordered)
// Simply give an empty map<string, string> mg (and a NULL pointer for tab_cpt) to skip including sequence to the transcript
GFFtranscript::GFFtranscript(vector<GFFentry> & vEntries, string sTAG, int shift, map<string, string> &mg, char *tab_cpt, long posStart, long posStop){
  this->init(vEntries, sTAG, shift, mg, tab_cpt, posStart, posStop);
}


void GFFtranscript::init(vector<GFFentry> & vEntries, string sTAG, int shift, map<string, string> &mg, char *tab_cpt, long posStart, long posStop){
  
  int i;

  sort(vEntries.begin(), vEntries.end(), my_GFFentry_compare);

  int nb = vEntries.size();
  string tSeq = "";

  this->posStart_ = -1;
  this->posStop_ = -1;

  this->id_ = vEntries[0].mat_[sTAG];
  this->source_ = vEntries[0].source_;

  this->chr_ = vEntries[0].seqID_;
  this->strand_ = vEntries[0].strand_;

  bool reverse = false;
  if( this->strand_=="-" ){
    reverse = true;
  }


  for(i=0 ; i<nb ; i++){

    GFFentry gff = vEntries[i];

    if( mg.size()>0 ){
      string stmp = mg[gff.seqID_].substr( (gff.start_ - shift) , (gff.end_ - gff.start_)+1 );
      if( reverse==true ){
	tSeq = stmp + tSeq;
      } else {
	tSeq = tSeq + stmp;
      }
    }

    pair<long, long> p;
    p.first = gff.start_ - shift;
    p.second = gff.end_ - shift;
    this->vexonsCoords_.push_back(p);
  }


  if(this->strand_=="-"){
    tSeq = my_reverse_cpt(tSeq, tab_cpt);
  }
  this->seq_ = tSeq;


  long start, end;
  this->getGenomicMinMaxPos(&start, &end);
  this->posMin_ = start;
  this->posMax_ = end;

}


// THIS FUNCTION IS NOT FINISHED -> TO DO !!!
string GFFtranscript::getIntron(int intronNumber){

  if(this->seq_.size()<1){
    fprintf(stderr, "WARNING: No sequence data, cannot get intron sequence!\n");
    return "";
  }

  int nbe = this->vexonsCoords_.size();

  if( intronNumber >= (nbe-1) ){
    fprintf(stderr, "WARNING: no intron (%d)\n", intronNumber);
    return "";
  }


}

void GFFtranscript::getAllEEJcoords(vector<long> &vCoords){

  bool cpt = false;
  if(this->strand_=="-"){ cpt = true; }
  long posEEJ = -1;

  int nbe = this->vexonsCoords_.size();

  for(int ee=1 ; ee<nbe ; ee++){
    if( cpt==false ){
      posEEJ = this->genomicToSplicedTranscript(this->vexonsCoords_[ee].first) - 1;
    } else {
      posEEJ = this->genomicToSplicedTranscript(this->vexonsCoords_[ee-1].first);
    }
    vCoords.push_back(posEEJ);
  }

}

void GFFtranscript::getAllIntronsCoordsInUnSplicedTranscripts(vector< pair<long, long> > &vCoords){

  long intronStart, intronEnd;
  bool cpt = false;

  if(this->strand_=="-"){ cpt = true; }

  int nbe = this->vexonsCoords_.size();

  for(int ee=1 ; ee<nbe ; ee++){

    pair<long,long> currentExon = this->vexonsCoords_[ee];
    pair<long,long> previousExon = this->vexonsCoords_[ee-1];

    if(cpt==false){
      intronStart = this->genomicToUnSplicedTranscript( previousExon.second ) + 1;
      intronEnd = this->genomicToUnSplicedTranscript( currentExon.first ) - 1;
    } else {
      intronStart = this->genomicToUnSplicedTranscript( previousExon.first ) + 1;
      intronEnd = this->genomicToUnSplicedTranscript( currentExon.second ) - 1;
    }

    pair<long,long> pp;
    pp.first = intronStart;
    pp.second = intronEnd;
    vCoords.push_back(pp);
  }

}


// THIS FUNCTION IS NOT FINISHED -> TO DO !!!
void GFFtranscript::getAllIntrons(vector<string> & vIntrons){

  if(this->seq_.size()<1){
    fprintf(stderr, "WARNING: No sequence data, cannot get intron sequence!\n");
    return;
  }

  int nbe = this->vexonsCoords_.size();

  for(int i=1 ; i<nbe ; i++){
    string intronSeq = getIntron(i-1);
    vIntrons.push_back(intronSeq);
  }

}


string GFFtranscript::getIntron(int intronNumber, map<string, string> &mg, char * tab_cpt){

  int nbe = this->vexonsCoords_.size();

  if( intronNumber >= (nbe-1) ){
    fprintf(stderr, "WARNING: no intron (%d)\n", intronNumber);
    return "";
  }

  long intronStart, intronEnd;

  bool reverse = true;
  if(this->strand_=="+"){ reverse = false; }

  if(reverse == false){
    intronStart = this->vexonsCoords_[intronNumber].second + 1;
    intronEnd = this->vexonsCoords_[intronNumber+1].first - 1;
  } else {
    intronStart = this->vexonsCoords_[intronNumber+1].second + 1;
    intronEnd = this->vexonsCoords_[intronNumber].first - 1;
  }

  long intronLg = (intronEnd - intronStart) + 1;
  string intronSeq = mg[this->chr_].substr(intronStart, intronLg);
  if(reverse==true){
    intronSeq = my_reverse_cpt(intronSeq, tab_cpt);
  }

  return intronSeq;

}


void GFFtranscript::getAllIntrons(vector<string> & vIntrons, map<string, string> &mg, char *tab_cpt){

  int nbe = this->vexonsCoords_.size();

  for(int i=1 ; i<nbe ; i++){
    string intronSeq = getIntron(i-1, mg, tab_cpt);
    vIntrons.push_back(intronSeq);
  }

}


string GFFtranscript::getExon(int exonNumber, map<string, string> &mg, char *tab_cpt){

  int nbe = this->vexonsCoords_.size();

  if(exonNumber >= nbe){
    fprintf(stderr, "WARNING: exon number too high (%d)\n", exonNumber);
    return "";
  }

  bool reverse = true;
  if(this->strand_=="+"){ reverse = false; }

  long exonStart = this->vexonsCoords_[exonNumber].first;
  long exonEnd = this->vexonsCoords_[exonNumber].second;

  long exonSize = (exonEnd - exonStart) + 1;

  string exonSeq = mg[this->chr_].substr(exonStart, exonSize);

  if(reverse==true){
    exonSeq = my_reverse_cpt(exonSeq, tab_cpt);
  }

  return exonSeq;

}


void GFFtranscript::getAllExons(vector<string> &vExons, map<string, string> &mg, char *tab_cpt){

  int nbe = this->vexonsCoords_.size();

  for(int i=0 ; i<nbe ; i++){
    string exonSeq = getExon(i, mg, tab_cpt);
    vExons.push_back(exonSeq);
  }

}


void GFFtranscript::markExonsGenomicScale(vector<bool> &v, bool mark){

  int nbe = this->vexonsCoords_.size();

  for(int i=0 ; i<nbe ; i++){
    long exonStart = this->vexonsCoords_[i].first;
    long exonEnd = this->vexonsCoords_[i].second;
    for(long ll=exonStart ; ll<=exonEnd ; ll++){
      v[ll] = mark;
    }
  }

}



void GFFtranscript::markIntronsGenomicScale(vector<bool> &v, bool mark){

  int nbe = this->vexonsCoords_.size();

  bool reverse = true;
  if(this->strand_=="+"){ reverse = false; }

  long intronStart, intronEnd;

  for(int i=1 ; i<nbe ; i++){

    if(reverse == false){
      intronStart = this->vexonsCoords_[i-1].second + 1;
      intronEnd = this->vexonsCoords_[i].first - 1;
    } else {
      intronStart = this->vexonsCoords_[i].second + 1;
      intronEnd = this->vexonsCoords_[i-1].first - 1;
    }

    for(long ll=intronStart ; ll<=intronEnd ; ll++){
      v[ll] = mark;
    }

  }

}

long GFFtranscript::getSplicedSize(void){

  long size = 0;

  int nbe = this->vexonsCoords_.size();
  for(int ee=0 ; ee<nbe ; ee++){
    pair<long,long> pp = this->vexonsCoords_[ee];
    long exonSize = (pp.second-pp.first) + 1;
    size += exonSize;
  }

  return size;

}

long GFFtranscript::splicedToGenomic_deprecated(long pos){

  int i;
  long sizePreviousExons = 0;
  int nbe = this->vexonsCoords_.size();

  long start, end;
  this->getGenomicMinMaxPos(&start, &end);

  long res, sizeExon;
  pair<long, long> pp;

  if( this->strand_ == "+" ){
    for(i=0 ; i<nbe ; i++){
      pp = this->vexonsCoords_[i];

      pp.first = this->genomicToSplicedTranscript(pp.first);
      pp.second = this->genomicToSplicedTranscript(pp.second);

      if(pos>=pp.first && pos<=pp.second){
	res = start + pos + sizePreviousExons;
	return res;
      }

      sizeExon = (pp.second - pp.first) + 1;
      sizePreviousExons += sizeExon;
    }
  }

  if( this->strand_ == "-" ){
    for(i=0 ; i<nbe ; i++){
      pp = this->vexonsCoords_[i];

      pp.first = this->genomicToSplicedTranscript(pp.first);
      pp.second = this->genomicToSplicedTranscript(pp.second);

      if(pos>=pp.second && pos<=pp.first){
	res = end - (pos + sizePreviousExons);
	return res;
      }

      sizeExon = (pp.second - pp.first) + 1;
      sizePreviousExons += sizeExon;
    }    
  }

  return -1;

}


long GFFtranscript::splicedToGenomic(long pos){

  int i;
  long sizePreviousExons = 0;
  long sizePreviousIntrons = 0;
  int nbe = this->vexonsCoords_.size();

  long start, end;
  this->getGenomicMinMaxPos(&start, &end);

  long res, sizeExon, sizeIntron;
  pair<long, long> pp;
  pair<long,long> pp_prev;

  if( this->strand_ == "+" ){
    for(i=0 ; i<nbe ; i++){
      pp = this->vexonsCoords_[i];

      sizeExon = (pp.second - pp.first) + 1;
      sizePreviousExons += sizeExon;

      if(i>0){
	pp_prev = this->vexonsCoords_[i-1];
	sizeIntron = (pp.first-1) - (pp_prev.second+1) + 1;
	sizePreviousIntrons += sizeIntron;
      }

      if(sizePreviousExons>pos){
	res = start + sizePreviousIntrons + pos;
	return res;
      }
      
    }
  }

  if( this->strand_ == "-" ){
    for(i=0 ; i<nbe ; i++){
      pp = this->vexonsCoords_[i];
      //fprintf(stderr, "%d - %d\n", pp.first, pp.second);

      sizeExon = (pp.second - pp.first) + 1;
      sizePreviousExons += sizeExon;

      if(i>0){
	pp_prev = this->vexonsCoords_[i-1];
	sizeIntron = (pp_prev.first-1) - (pp.second+1) + 1;
	sizePreviousIntrons += sizeIntron;
      }

      if(sizePreviousExons>pos){
	res = end - (pos+sizePreviousIntrons);
	return res;
      }

    }
  }

  return -1;

}


long GFFtranscript::genomicToSplicedTranscript(long pos){

  int i;
  long sizePreviousExons = 0;
  int nbe = this->vexonsCoords_.size();

  long res, sizeExon;
  pair<long, long> pp;

  if( this->strand_ == "+" ){
    for(i=0 ; i<nbe ; i++){
      pp = this->vexonsCoords_[i];
      if(pos>=pp.first && pos<=pp.second){
	res = (pos - pp.first) + sizePreviousExons;
	return res;
      }

      sizeExon = (pp.second - pp.first) + 1;
      sizePreviousExons += sizeExon;
    }
  }

  if( this->strand_ == "-" ){
    for(i=0 ; i<nbe ; i++){
      pp = this->vexonsCoords_[i];
      if(pos>=pp.first && pos<=pp.second){
	res = (pp.second - pos) + sizePreviousExons;
	return res;
      }

      sizeExon = (pp.second - pp.first) + 1;
      sizePreviousExons += sizeExon;
    }
    
  }
  
  return -1;
}

int GFFtranscript::allSplicedToGenomic(map<int,int> &m){

  int i, j;
  long sizePreviousExons = 0;
  int nbe = this->vexonsCoords_.size();

  long res, sizeExon;
  pair<long, long> pp;

  if( this->strand_ == "+" ){
    for(i=0 ; i<nbe ; i++){
      pp = this->vexonsCoords_[i];
      sizeExon = (pp.second - pp.first) + 1;
      for(j=0 ; j<sizeExon ; j++){
	m[j+sizePreviousExons] = pp.first + j;
      }      
      sizePreviousExons += sizeExon;
    }
  }

  if( this->strand_ == "-" ){
    for(i=0 ; i<nbe ; i++){
      pp = this->vexonsCoords_[i];
      sizeExon = (pp.second - pp.first) + 1;
      for(j=0 ; j<sizeExon ; j++){
	m[j+sizePreviousExons] = pp.second - j;
      }
      sizePreviousExons += sizeExon;
    }
  }
 
  return 0;  
}


long GFFtranscript::genomicToUnSplicedTranscript(long pos){

  int i;
  long sizePreviousExons = 0;
  int nbe = this->vexonsCoords_.size();

  long res, sizeExon;
  pair<long, long> ppCurrentExon;
  pair<long, long> ppPreviousExon;

  if( this->strand_ == "+" ){
    for(i=0 ; i<nbe ; i++){

      ppCurrentExon = this->vexonsCoords_[i];
      ppPreviousExon = this->vexonsCoords_[i-1];

      if(i>0){
	int sizePreviousIntron = ( (ppCurrentExon.first-1) - (ppPreviousExon.second+1) ) + 1;
	sizePreviousExons += sizePreviousIntron;
      }


      if(pos>=ppCurrentExon.first && pos<=ppCurrentExon.second){
	res = (pos - ppCurrentExon.first) + sizePreviousExons;
	return res;
      }

      sizeExon = (ppCurrentExon.second - ppCurrentExon.first) + 1;
      sizePreviousExons += sizeExon;
    }
  }

  if( this->strand_ == "-" ){
    for(i=0 ; i<nbe ; i++){

      ppCurrentExon = this->vexonsCoords_[i];
      ppPreviousExon = this->vexonsCoords_[i-1];

      if(i>0){
	int sizePreviousIntron = ( (ppPreviousExon.first-1) - (ppCurrentExon.second+1) ) + 1;
	sizePreviousExons += sizePreviousIntron;
      }


      if(pos>=ppCurrentExon.first && pos<=ppCurrentExon.second){
	res = (ppCurrentExon.second - pos) + sizePreviousExons;
	return res;
      }

      sizeExon = (ppCurrentExon.second - ppCurrentExon.first) + 1;
      sizePreviousExons += sizeExon;
    }
    
  }

  
  return -1;

}

void GFFtranscript::getSpliceUnSpliceAndFlanking(int intronNumber, int flankingSize, map<string, string> &mg, char *tab_cpt, string &spliced, string &unSpliced){

  long start, end, lg;
  string flank5p, flank3p;

  string intronSeq = this->getIntron(intronNumber, mg, tab_cpt);
  this->getFlankingFromIntron(intronNumber, flankingSize, mg, tab_cpt, flank5p, flank3p);

  unSpliced = flank5p + intronSeq + flank3p;
  spliced = flank5p + flank3p;

}


void GFFtranscript::getFlankingFromIntron(int intronNumber, int flankingSize, map<string, string> &mg, char *tab_cpt, string &flank5p, string &flank3p){

  long start, end, lg;
  //string flank5p, flank3p;

  // First, I extract the 5' flanking exonic sequence
  if(this->strand_ == "+"){
    end = this->genomicToSplicedTranscript( this->vexonsCoords_[intronNumber].second );
  } else {
    end = this->genomicToSplicedTranscript( this->vexonsCoords_[intronNumber].first );
  }

  start = (end - flankingSize) + 1;
  if(start<0){ start = 0; }

  lg = (end - start) + 1;
  flank5p = this->seq_.substr(start, lg);

  if(this->strand_ == "+"){
    start = this->genomicToSplicedTranscript( this->vexonsCoords_[intronNumber+1].first );
  } else {
    start = this->genomicToSplicedTranscript( this->vexonsCoords_[intronNumber+1].second );
  }

  end = (start + flankingSize) - 1;
  if( end>this->seq_.size() ){ end = this->seq_.size(); }

  lg = (end - start) + 1;
  flank3p = this->seq_.substr(start, lg);

}



string GFFtranscript::getIntronWithFlanking(int intronNumber, int flankingSize, map<string, string> &mg, char *tab_cpt){

  string flank5p, flank3p;

  string intronSeq = this->getIntron(intronNumber, mg, tab_cpt);
  this->getFlankingFromIntron(intronNumber, flankingSize, mg, tab_cpt, flank5p, flank3p);

  string s = flank5p + intronSeq + flank3p;
  return s;
}



bool GFFtranscript::intronIsNeverExonic(int intronNumber, map<string, vector<bool> > &mpi, bool exonic){


  map<string, vector<bool> >::iterator it = mpi.find(this->chr_);

  long start, end;

  if(this->strand_ == "+"){
    start = this->vexonsCoords_[intronNumber].second + 1;
    end = this->vexonsCoords_[intronNumber+1].first - 1;
  } else {
    start = this->vexonsCoords_[intronNumber+1].second + 1;
    end = this->vexonsCoords_[intronNumber].first - 1;
  }

  for(long ll=start ; ll<=end ; ll++){
    if( (it->second).at(ll) == exonic ){
      return false;
    }
  }

  return true;
}


int GFFtranscript::getNbExons(void){
  return this->vexonsCoords_.size();
}


int GFFtranscript::getNbIntrons(void){
  return (this->vexonsCoords_.size() - 1);
}



void GFFtranscript::printCoords(FILE *nf){

  fprintf(nf, "%s\t%s\t%s", 
	  this->chr_.c_str(),
	  this->id_.c_str(),
	  this->strand_.c_str()
	  );
  
  vector< pair<long, long> > ve = this->vexonsCoords_;
  int nbe = ve.size();
  for(int i=0 ; i<nbe ; i++){
    fprintf(nf, "\t[%ld - %ld]", ve[i].first, ve[i].second);
  }
  
  fprintf(nf, "\n");  

}



string GFFtranscript::getUnSplicedSeq(map<string, string> &mg, char *tab_cpt){

  string exon, intron, stmp;
  string seq = "";

  int nbi = this->vexonsCoords_.size() - 1;

  int i;
  for(i=0 ; i<nbi ; i++){
    exon = this->getExon(i, mg, tab_cpt);
    intron = this->getIntron(i, mg, tab_cpt);
    stmp = exon + intron;

    seq = seq + stmp;

    // getExon

    // getIntron

    // concatenate both sequences

  }

  exon = this->getExon(i, mg, tab_cpt);
  seq = seq + exon;

  return seq;
}


string GFFtranscript::getSplicedSeq(void){
  return this->seq_;
}


void GFFtranscript::getGenomicMinMaxPos(long *start, long *end){

  int nbExons = this->vexonsCoords_.size();

  if( this->strand_ == "+" ){
    (*start) = this->vexonsCoords_[0].first;
    (*end) = this->vexonsCoords_[nbExons-1].second;
  } else {
    (*start) = this->vexonsCoords_[nbExons-1].first;
    (*end) = this->vexonsCoords_[0].second;
  }

}

string GFFtranscript::getUpStreamGenomicSeq(int nbBases, string &chrSeq, char *tab_cpt, int *nbDone){

  long start, end;
  this->getGenomicMinMaxPos(&start, &end);

  string seq;
  int chrStart, lg;

  if(this->strand_=="+"){
    chrStart = start - nbBases;
    if(chrStart<0){
      chrStart = 0;
    }

    lg = start - chrStart;
    seq = chrSeq.substr(chrStart, lg);
    (*nbDone) = lg;

  } else {

    chrStart = end+1;
    lg = nbBases;
    int margin = chrSeq.size() - chrStart;
    if( margin<nbBases ){
      lg = margin;
    }

    seq = chrSeq.substr(chrStart, lg);
    seq = my_reverse_cpt(seq, tab_cpt);
    (*nbDone) = lg;

  }


  return seq;
}



string GFFtranscript::getDownStreamGenomicSeq(int nbBases, string &chrSeq, char *tab_cpt, int *nbDone){

  long start, end;
  this->getGenomicMinMaxPos(&start, &end);

  string seq;
  int chrStart, lg;

  if(this->strand_=="+"){

    chrStart = end + 1;

    int margin = chrSeq.size() - chrStart;
    if( margin<nbBases ){
      lg = margin - 1;
    } else {
      lg = nbBases;
    }

    seq = chrSeq.substr(chrStart, lg);
    (*nbDone) = lg;

  } else {

    chrStart = start - nbBases;

    if(chrStart<0){
      chrStart = 0;
    }

    lg = start - chrStart;
    seq = chrSeq.substr(chrStart, lg);
    (*nbDone) = lg;

    seq = my_reverse_cpt(seq, tab_cpt);

  }

  return seq;


}

bool GFFtranscript::isOnPlusStrand(void){
  bool plusStrand = true;
  if(this->strand_=="-"){
    plusStrand = false;
  }
  return plusStrand;
}



bool my_GFFentry_compare(GFFentry e1, GFFentry e2){

  if(e1.strand_=="+"){
    if( e1.start_ < e2.start_ ){
      return true;
    } else {
      return false;
    }
  } else {
    if( e1.start_ > e2.start_ ){
      return true;
    } else {
      return false;
    }
  }

  return false;
}

