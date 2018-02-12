/*

indelCoverage.cxx

 */


#include "indelCoverage.hxx"

using namespace std;

C_indelCoverage::C_indelCoverage(void){
  this->cov_ = 0;
  this->ins_ = 0;
  this->del_ = 0;
}


C_indelCoverage::C_indelCoverage(long cov, long ins, long del){

  this->cov_ = cov;
  this->ins_ = ins;
  this->del_ = del;

}


void C_indelCoverage::update(long cov, long ins, long del){
  this->cov_ += cov;
  this->ins_ += ins;
  this->del_ += del;
}


void C_indelCoverage::display(FILE *nf){

  fprintf(nf, 
	  "%d\t%d\t%d",
	  this->cov_,
	  this->ins_,
	  this->del_
	  );

}
