/*

parse_blast.cxx

Some utility functions to parse a BLAST output file.

 */

#ifndef _PARSE_BLAST_CXX_
#define _PARSE_BLAST_CXX_

#include <stdlib.h>

#include "parse_blast.hxx"


using namespace std;


blast_m8::blast_m8(void){
  this->queryId_ = "";
  this->subjectId_ = "";
}


blast_m8::blast_m8(string &s){

  this->make(s);

}

void blast_m8::make(string &s){

  string stmp;
  istringstream iss(s);

  getline(iss, this->queryId_, '\t');
  getline(iss, this->subjectId_, '\t');

  getline(iss, stmp, '\t');
  this->prctIdent_ = atof(stmp.c_str());

  getline(iss, stmp, '\t');
  this->alnLength_ = atoi(stmp.c_str());

  getline(iss, stmp, '\t');
  this->mismatchCount_ = atoi(stmp.c_str());

  getline(iss, stmp, '\t');
  this->gapOpenCount_ = atoi(stmp.c_str());

  getline(iss, stmp, '\t');
  this->queryStart_ = atoi(stmp.c_str());

  getline(iss, stmp, '\t');
  this->queryEnd_ = atoi(stmp.c_str());

  getline(iss, stmp, '\t');
  this->subjectStart_ = atoi(stmp.c_str());

  getline(iss, stmp, '\t');
  this->subjectEnd_ = atoi(stmp.c_str());

  getline(iss, stmp, '\t');
  this->strEval_ = stmp;
  this->eVal_ = atof(stmp.c_str());

  getline(iss, stmp, '\t');
  this->bitScore_ = atof(stmp.c_str());

}


void blast_m8::print(FILE *nf){

  fprintf(
	  nf, 
	  "%s\t%s\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%.1f\n",
	  this->queryId_.c_str(),
	  this->subjectId_.c_str(),
	  this->prctIdent_,
	  this->alnLength_,
	  this->mismatchCount_,
	  this->gapOpenCount_,
	  this->queryStart_,
	  this->queryEnd_,
	  this->subjectStart_,
	  this->subjectEnd_,
	  this->strEval_.c_str(),
	  this->bitScore_
	  );

}

#endif
