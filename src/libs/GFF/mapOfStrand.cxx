#ifndef _mapOfStrand_CXX_
#define _mapOfStrand_CXX_

/*

mapOfStrand.cxx

 */


#include "mapOfStrand.hxx"



mapOfStrand::mapOfStrand(void){
  ;
}


mapOfStrand::mapOfStrand(map<string, string> &mg, map<string, GFFtranscript> &mt){

  this->init(mg);

  map<string, GFFtranscript>::iterator it;
  for(it=mt.begin() ; it!=mt.end() ; it++){
    long min, max;
    GFFtranscript gt = it->second;
    gt.getGenomicMinMaxPos(&min, &max);
    string chrID = gt.chr_;

    char cs = '.';
    if(gt.strand_=="+"){
      cs = '+';
    } else {
      cs = '-';
    }

    //fprintf(stdout, "%s\t%d\t%d\n", gt.strand_.c_str(), min, max);

    for(long pos=min ; pos<=max ; pos++){
      char currentPos = this->m_[chrID].at(pos);
      if( currentPos=='.' ){
	this->m_[chrID][pos] = cs;
      } else {
	if( currentPos != cs ){
	  this->m_[chrID][pos] = 'b';
	}
      }
    }

  }

}


void mapOfStrand::init(map<string, string> &mg){

  map<string, string>::iterator it;
  for(it=mg.begin() ; it!=mg.end() ; it++){
    string chrID = it->first;
    long lg = (it->second).size();
    //vector<char> vtmp('.', lg);
    //long nb = vtmp.size();
    for(long ll=0 ; ll<lg ; ll++){
      this->m_[chrID].push_back('.');
    }
  }

}



char mapOfStrand::getStrand(string chrID, long pos){

  map<string, vector<char> >::iterator it = this->m_.find(chrID);

  if( it==this->m_.end() ){
    return '?';
  }

  long vs = it->second.size();
  if( pos<0 || pos>=vs )
    return '?';

  return it->second.at(pos);
  //(*c) = it->at(pos);
  //return true;
  //return this->m_[chrID].at(pos);
}


#endif
