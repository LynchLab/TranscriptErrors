#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>


#include "naive_utils.hxx"

using namespace std;


bool fillVectorFromArgs(vector<string> &vArgs, string &sArgs, char sep){

  string stmp;

  istringstream iss(sArgs);
  while( getline(iss, stmp, sep) ){
    vArgs.push_back(stmp);
  }

  return true;

}


int utils_findPosInMate(ifstream &ficBlatPos, string &id){

  string line, local_id, stmp;
  int csLg, csStart;

  while( getline(ficBlatPos, line) ){
    istringstream iss(line);
    getline(iss, local_id, '\t');

    if( local_id == id ){
      getline(iss, stmp, '\t');
      getline(iss, stmp, '\t');
      csStart = atoi(stmp.c_str());
      return csStart;
    }
  }

  return -1;
}


int naive_utils_getDebreakedCs(
			       CS_consensus &cs, 
			       string &id, 
			       ifstream &ficBlatPos,
			       ifstream &ficCsPos,
			       ifstream &ficReads, 
			       ifstream &ficMates, 
			       bool paired, 
			       char *tab_cpt
			       ){

  string lineBlat;
  int csLgParam, alnLgParam, breakPosParam;
  return naive_utils_getDebreakedCs(cs, id, ficBlatPos, lineBlat, ficCsPos, ficReads, ficMates, paired, tab_cpt, &csLgParam, &alnLgParam, &breakPosParam);

}



int naive_utils_getDebreakedCs(
			       CS_consensus &cs, 
			       string &id, 
			       ifstream &ficBlatPos, string &lineBlat,
			       ifstream &ficCsPos,
			       ifstream &ficReads, 
			       ifstream &ficMates, 
			       bool paired, 
			       char *tab_cpt,
			       int *csLgParam,
			       int *alnLgParam,
			       int *breakPosParam
			       ){

  string stmp, blatId;
  int csLg, alnLg, breakPos;

  blatId= "";

  while( blatId!=id && getline(ficBlatPos, lineBlat) ){
    istringstream iss(lineBlat);
    getline(iss, blatId, '\t');

    if(blatId==id){
      getline(iss, stmp, '\t');
      csLg = atoi(stmp.c_str());
      getline(iss, stmp, '\t');
      alnLg = atoi(stmp.c_str());
      getline(iss, stmp, '\t');
      breakPos = atoi(stmp.c_str());

      (*csLgParam) = csLg;
      (*alnLgParam) = alnLg;
      (*breakPosParam) = breakPos;

      // DEBOG FOR STRANGE CASES THAT CAUSE PROBLEM
      if( alnLg>=csLg && breakPos>csLg ){
	return 1;
      }

    }

  }

  if( blatId!=id ){
    fprintf(stderr, "ERROR: '%s' not found among blat hits.\n", id.c_str());
    return -1;
  }


  fastqRead read, mate;

  int posCsInMate;
  if( paired==true ){
    //fprintf(stderr, "Searching for '%s' in '%s'\n", id.c_str(), argv[ARG_FIC_BLAT_
    posCsInMate = utils_findPosInMate(ficCsPos, id);
    if( posCsInMate<0 ){
      fprintf(stderr, "ERROR: %s not found in file of positions.\n", id.c_str());
      return -2;
    }
  }

  while( read.readNext(ficReads) && read.id_ != id ){
    if( paired==true ){
      mate.readNext(ficMates);
    }
  }    
  if( paired==true ){
    mate.readNext(ficMates);
  }


  if( read.id_ != id ){
    cerr<<"ERROR: "
	<<id
	<<" was not found in the reads file."
	<<endl;
    return -3;
  }

  if( paired==true && read.id_ != mate.id_){
    cerr<<"ERROR: read and mate id do not match ('"
	<<read.id_<<"' and '"<<mate.id_<<"')"
	<<endl;
    return -4;
  }


  read.trimLowQuality('#');
  if( paired==true ){
    mate.trimLowQuality('#');
    string mateSeq = my_reverse_cpt(mate.seq_, tab_cpt);
    string mateQual = my_reverse(mate.squal_);
    cs.make(read.seq_, read.squal_, mateSeq, mateQual, csLg, breakPos, posCsInMate);
    //cs.make(read.seq_, read.squal_, mateSeq, mateQual, debreakedRead.seq_, minID);
  } else {
    ;
    // TO DO !!!
    cs.make(read.seq_, read.squal_, csLg, breakPos);
    //cs.make(read.seq_, read.squal_, breakPos, minID);
    //cs.make(read.seq_, read.squal_, debreakedRead.seq_, minID);
  }

  //fprintf(stdout, ">%s\n", read.id_.c_str());
  //cs.prettyPrint(stdout, false);

  return 0;
  
}

