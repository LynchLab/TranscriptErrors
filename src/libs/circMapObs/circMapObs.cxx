/*

circMapObs.cxx

 */


#include <stdlib.h>

#include <iostream>
#include <map>

#include "circMapObs.hxx"


using namespace std;



circMapObs::circMapObs(void){
  ;
}


circMapObs::circMapObs(std::map<std::string,std::string> &mg){
  this->init(mg);
}

circMapObs::circMapObs(std::vector<string> & vFiles){
  this->init(vFiles);
}

circMapObs::circMapObs(char *fname){
  this->init(fname);
}


void circMapObs::init(std::vector<string> & vFiles){

  int nb = vFiles.size();

  if(nb==0)
    return;

  this->init((char *)vFiles[0].c_str());

  for(int i=1 ; i<nb ; i++){
    this->add((char *)vFiles[i].c_str());
  }

}


void circMapObs::init(std::map<std::string,std::string> &mg){
  map<string,string>::iterator itm;
  for(itm=mg.begin() ; itm!=mg.end() ; itm++){
    string chrID = itm->first;
    int chrLg = (itm->second).size();
    //vector<long> vtmp(chrLg, 0);
    //this->mObs_[chrID] = vtmp;
    //map<string, vector<long> >::iterator itt = this->mObs_.find(chrID);
    for(int i=0 ; i<chrLg ; i++){
      this->mObs_[chrID].push_back(0);
      //(itt->second).push_back(0);
    }
  }
}


void circMapObs::init(char *fname){
  string line, chrID, sPos, sObs;
  long nbObs;

  map<string, vector<long> >::iterator it;
  string prev_chrID = "";

  this->mObs_.clear();

  ifstream fic(fname);
  while( getline(fic, line) ){
    istringstream iss(line);

    getline(iss, chrID, '\t');
    getline(iss, sPos, '\t');
    getline(iss, sObs, '\t');


    nbObs = atoi(sObs.c_str());
    this->mObs_[chrID].push_back(nbObs);
  }
}


void circMapObs::add(char *fname){
  string line, chrID, sPos, sObs;
  int pos;
  long nbObs;

  //this->mObs_.clear();

  ifstream fic(fname);
  while( getline(fic, line) ){
    istringstream iss(line);

    getline(iss, chrID, '\t');
    getline(iss, sPos, '\t');
    getline(iss, sObs, '\t');

    int pos = atoi(sPos.c_str());
    nbObs = atoi(sObs.c_str());

    this->mObs_[chrID].at(pos) += nbObs;
  }  
}


void circMapObs::addObs(std::string chrID, int pos){
  (this->mObs_[chrID].at(pos))++;
}

void circMapObs::addObs(std::string chrID, vector<long> vpos){
  map<string, vector<long> >::iterator it = this->mObs_.find(chrID);

  vector<long>::iterator itv;
  for(itv=vpos.begin() ; itv!=vpos.end() ; itv++){
    (it->second).at((*itv))++;
  }

  //(this->mObs_[chrID].at(pos))++;
}


void circMapObs::setObs(std::string chrID, int pos, long nb){
  this->mObs_[chrID].at(pos) = nb;
}

long circMapObs::getObs(std::string chrID, int pos){
  return this->mObs_[chrID].at(pos);
}


long circMapObs::getSum(void){

  long sum = 0;

  map<string, vector<long> >::iterator itm;
  for(itm=this->mObs_.begin() ; itm!=this->mObs_.end() ; itm++){
    vector<long>::iterator itv;
    for(itv=(itm->second).begin() ; itv!=(itm->second).end() ; itv++){
      sum += (*itv);
    }
  }

  return sum;

}


void circMapObs::print(FILE *nf){
  map<string,vector<long> >::iterator itm;
  for(itm=this->mObs_.begin() ; itm!=this->mObs_.end() ; itm++){
    string chrID = itm->first;
    int chrLg = (itm->second).size();
    for(int i=0 ; i<chrLg ; i++){
      fprintf(nf, "%s\t%d\t%ld\n", chrID.c_str(), i, this->mObs_[chrID].at(i));
    }
  }
}



void circMapObs::printWithBase(FILE *nf, std::map<std::string,std::string> & mg, mapOfStrand &mos){

  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);

  map<string,vector<long> >::iterator itm;
  for(itm=this->mObs_.begin() ; itm!=this->mObs_.end() ; itm++){
    string chrID = itm->first;
    int chrLg = (itm->second).size();

    for(int i=0 ; i<chrLg ; i++){

      char genomicBase = mg[chrID].at(i);
      char cStrand = mos.getStrand(chrID, i);

      if(cStrand!='+' && cStrand!='-'){
	genomicBase = '?';
      } else {
	if(cStrand=='-'){
	  genomicBase = tab_cpt[genomicBase];
	}
      }
      fprintf(nf, "%s\t%d\t%ld\t%c\n", chrID.c_str(), i, this->mObs_[chrID].at(i), genomicBase);
    }
  }
}


circMapObsAndBase::circMapObsAndBase(void){
  ;
}


circMapObsAndBase::circMapObsAndBase(char *fname){
  this->init(fname);
}


void circMapObsAndBase::init(char *fname){
  string line, chrID, stmp;

  ifstream fic(fname);
  while( getline(fic, line) ){
    istringstream iss(line);
    getline(iss, chrID, '\t');
    getline(iss, stmp, '\t');
    long gPos = atol(stmp.c_str());

    getline(iss, stmp, '\t');
    long nbObs = atol(stmp.c_str());

    getline(iss, stmp, '\t');
    char base = stmp[0];

    pair<long,char> pp(nbObs, base);
    this->m_[chrID][gPos] = pp;

  }

}

int circMapObsAndBase::get(std::string chrID, long gPos, std::pair<long,char> &pp){

  map<string, map<long, pair<long,char> > >::iterator itm = this->m_.find(chrID);
  if( itm==this->m_.end() ){
    return -1;
  }

  map<long, pair<long,char> >::iterator itp = itm->second.find(gPos);
  if( itp==itm->second.end() ){
    return -2;
  }

  pp.first = itp->second.first;
  pp.second = itp->second.second;

  return 0;
}
