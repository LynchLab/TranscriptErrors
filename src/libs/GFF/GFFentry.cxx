/*

GFFentry.cxx

 */


#include <stdlib.h>

#include "GFFentry.hxx"


GFFentry::GFFentry(void){
  ;
}



GFFentry::GFFentry(string line, string format){

  string stmp, sattributes, scouple, stag, svalue;

  istringstream iss(line);

  getline(iss, this->seqID_, '\t');
  getline(iss, this->source_, '\t');
  getline(iss, this->feature_, '\t');

  getline(iss, stmp, '\t');
  this->start_ = atol(stmp.c_str());
  getline(iss, stmp, '\t');
  this->end_ = atol(stmp.c_str());

  getline(iss, this->s_score_, '\t');
  getline(iss, this->strand_, '\t');

  getline(iss, stmp, '\t');
  this->phase_ = atoi(stmp.c_str());

  getline(iss, sattributes, '\t');
  istringstream issa(sattributes);

  while(getline(issa, scouple, ';')){
    istringstream issc(scouple);

    if(format=="GTF"){

      getline(issc, stag, ' ');
      if(stag==""){
	getline(issc, stag, ' ');
      }

      getline(issc, svalue);
      if(svalue[0]=='"'){
	svalue = svalue.substr(1, svalue.size()-2);
      }

    } else {
      getline(issc, stag, '=');
      getline(issc, svalue, '=');
      while(stag.size()>0 && stag[0]==' '){
	stag.erase(0, 1);
      }
    }

    if(stag.size()<1){
      fprintf(stderr, "WARNING: stag skipped because of size zero (line: '%s')\n", line.c_str());
    } else {
      this->mat_.insert(map<string,string>::value_type(stag, svalue));
    }

  }

}
