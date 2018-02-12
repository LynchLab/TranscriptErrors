#ifndef _CIRC_MAP_COV_SIMPLE_CXX_
#define _CIRC_MAP_COV_SIMPLE_CXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>

#include <stdlib.h>

#include "circMapIndels.hxx"


using namespace std;

C_indelPos::C_indelPos(void){
  this->cov_ = 0;
  this->nbIns_ = 0;
  this->nbDel_ = 0;
}


C_indelPos::C_indelPos(const long cov, const long nbIns, const long nbDel){
  this->cov_ = cov;
  this->nbIns_ = nbIns;
  this->nbDel_ = nbDel;
}


void C_indelPos::get(long *cov, long *nbIns, long *nbDel){
  (*cov) = this->cov_;
  (*nbIns) = this->nbIns_;
  (*nbDel) = this->nbDel_;
}

void C_indelPos::set(const long cov, const long nbIns, const long nbDel){
  this->cov_ = cov;
  this->nbIns_ = nbIns;
  this->nbDel_ = nbDel;
}

void C_indelPos::add(C_indelPos &ip){
  long cov = 0;
  long nbIns = 0;
  long nbDel = 0;
  ip.get(&cov, &nbIns, &nbDel);
  this->add(cov, nbIns, nbDel);
}

void C_indelPos::add(long cov, long nbIns, long nbDel){
  this->cov_ += cov;
  this->nbIns_ += nbIns;
  this->nbDel_ += nbDel;
}


void C_indelPos::addCov(void){
  this->cov_++;
}

void C_indelPos::addCov(const long nb){
  this->cov_ += nb;
}



void C_indelPos::addIns(void){
  this->nbIns_++;
}

void C_indelPos::addIns(const long nb){
  this->nbIns_ += nb;
}


void C_indelPos::addDel(void){
  this->nbDel_++;
}

void C_indelPos::addDel(const long nb){
  this->nbDel_ += nb;
}

void C_indelPos::print(FILE *nf){
  fprintf(
	  nf, 
	  "%ld\t%ld\t%ld", 
	  this->cov_,
	  this->nbIns_,
	  this->nbDel_
	  );
}


bool C_indelPos::isSuspect(const long minCoverage, const double maxFracAlter){

  long cov = this->cov_;
  if(cov<minCoverage)
    return true;

  double r = (double)((double)(this->nbIns_+this->nbDel_)/(double)cov);
  if( r>maxFracAlter ){
    return true;
  }
  return false;
  
}


#endif
