/*

demulti.cxx

Creates all possible combinations of a consensus sequences.

 */


#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <algorithm>

#include "UTILS_FASTA/utils_fasta.hxx"


using namespace std;


#define ARG_FILE_CS 1
#define ARG_FILE_IDS 2


void demult(fastqRead &fr, FILE *nf);

int main(int argc, char *argv[]){

  string line;

  vector<string> vids;
  ifstream fIDs(argv[ARG_FILE_IDS]);
  while( getline(fIDs, line) ){
    vids.push_back(line);
  }

  fastqRead fr;
  ifstream ficCs(argv[ARG_FILE_CS]);

  while( fr.readNext(ficCs) ){
    if( find(vids.begin(), vids.end(), fr.id_)!=vids.end() ){
      demult(fr, stdout);
    }
  }

  return 0;
}

void demult(fastqRead &fr, FILE *nf){

  int lg = fr.seq_.size();

  for(int i=0 ; i<lg ; i++){
    string seqTmp = fr.seq_.substr(i, lg-i) + fr.seq_.substr(0, i);
    string qualTmp = fr.squal_.substr(i, lg-i) + fr.squal_.substr(0, i);
    fprintf(nf, 
	    "@%s_V%d\n%s\n%c\n%s\n",
	    fr.fullID_.c_str(),
	    i,
	    seqTmp.c_str(),
	    '+',
	    qualTmp.c_str()
	    );
  }


}

