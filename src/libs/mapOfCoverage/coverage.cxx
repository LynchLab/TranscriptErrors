/*

coverage.cxx

 */

#ifndef _coverage_CXX_
#define _coverage_CXX_

#include "coverage.hxx"


coverage::coverage(void){
  ;
}


coverage::coverage(long val){
  for(int i=0 ; i<NB_BASES_MOC ; i++){
    this->adjustedCov_[i] = val;
  }
}


coverage::coverage(long *tabVal){
  for(int i=0 ; i<NB_BASES_MOC ; i++){
    this->adjustedCov_[i] = tabVal[i];
  }
}


void coverage::display(FILE* nf){
  for(int ii=0 ; ii<NB_BASES_MOC ; ii++){
    fprintf(nf, "\t%d", this->adjustedCov_[ii]);
  }
}


void coverage::add(long *tabVal){
  for(int ii=0 ; ii<NB_BASES_MOC ; ii++){
    this->adjustedCov_[ii] += tabVal[ii];
  }
}

void coverage::add(char c, char *tabBases){
  for(int ii=0 ; ii<NB_BASES_MOC ; ii++){
    if( tabBases[ii] == c ){
      (this->adjustedCov_[ii])++;
      return;
    }
  }
}


void coverage::addFromGenomic(coverage &cv, char *tabBases, char strand, char *tabCpt){

  int ii, jj;

  if( strand == _STRAND_PLUS_ ){
    for(ii=0 ; ii<NB_BASES_MOC ; ii++){
      this->adjustedCov_[ii] += cv.adjustedCov_[ii];
    }
  }


  if( strand == _STRAND_MINUS_ ){
    for(ii=0 ; ii<NB_BASES_MOC ; ii++){
      for(jj=0 ; jj<NB_BASES_MOC && tabBases[jj]!=tabCpt[tabBases[ii]] ; jj++)
	;
      if( jj < NB_BASES_MOC ){
	this->adjustedCov_[ii] += cv.adjustedCov_[jj];
      }
    }
  }

}


#endif
