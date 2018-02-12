/*

utils_fasta.cxx

Just a small library with some usefull function to handle fasta files.

 */

#ifndef _UTILS_FASTA_CXX_
#define _UTILS_FASTA_CXX_

#include "utils_fasta.hxx"

using namespace std;

void loadFasta(map<string, string> & m, char *fname, bool extractID, bool verbose, bool debug){

  if(verbose==true){
    fprintf(stderr, "Loading %s ... ", fname);
    fflush(stderr);
  }

  string ligne, sid, seq;
  bool preums = true;
  vector<string> v;

  ifstream fic(fname);
  seq = "";
  while( getline(fic, ligne) ){
    if(ligne.size()>1 && ligne[0]=='>'){
      if(preums==false){
	m.insert(map<string, string>::value_type(sid, seq));
	seq.erase();
      }

      int pos_pipe = ligne.find_first_of("|");
      int pos_blanc = ligne.find_first_of(" ");
      int pos_tab = ligne.find_first_of("\t");

      int pos_sep = pos_pipe;
      if(pos_blanc!=string::npos && (pos_sep==string::npos ||  pos_blanc<pos_sep)){ pos_sep = pos_blanc; }
      if(pos_tab!=string::npos && (pos_sep==string::npos || pos_tab<pos_sep)){ pos_sep = pos_tab; }

      if(extractID==true){
	//sid = ligne.substr(1, ligne.find(" ")-1);
	sid = ligne.substr(1, pos_sep-1);
      } else {
	sid = ligne.substr(1, ligne.size());
      }

      if(debug==true){
	fprintf(stderr, "'%s'\n", sid.c_str());
	fflush(stderr);
	//cout<<"'"<<sid<<"'"<<endl;
      }
    } else {
      //seq = seq + ligne;
      seq.append(ligne);
    }
    preums = false;
  }

  m.insert(map<string, string>::value_type(sid, seq));

  if(verbose==true){
    fprintf(stderr, "done.\n");
    fflush(stderr);
  }
  
}

void loadFasta(map<string, string> & m, char *fname, string chrID, bool extractID, bool verbose, bool debug){

  if(verbose==true){
    fprintf(stderr, "Loading %s ... ", fname);
    fflush(stderr);
  }

  string ligne, sid, seq;
  bool preums = true;
  vector<string> v;

  ifstream fic(fname);
  seq = "";
  while( getline(fic, ligne) ){
    if(ligne.size()>1 && ligne[0]=='>'){
      if(preums==false){

	if(sid==chrID){
	  m.insert(map<string, string>::value_type(sid, seq));
	}

	seq.erase();
      }

      int pos_pipe = ligne.find_first_of("|");
      int pos_blanc = ligne.find_first_of(" ");
      int pos_tab = ligne.find_first_of("\t");

      int pos_sep = pos_pipe;
      if(pos_blanc!=string::npos && (pos_sep==string::npos ||  pos_blanc<pos_sep)){ pos_sep = pos_blanc; }
      if(pos_tab!=string::npos && (pos_sep==string::npos || pos_tab<pos_sep)){ pos_sep = pos_tab; }

      if(extractID==true){
	//sid = ligne.substr(1, ligne.find(" ")-1);
	sid = ligne.substr(1, pos_sep-1);
      } else {
	sid = ligne.substr(1, ligne.size());
      }

      if(debug==true){
	fprintf(stderr, "'%s'\n", sid.c_str());
	fflush(stderr);
	//cout<<"'"<<sid<<"'"<<endl;
      }
    } else {
      //seq = seq + ligne;
      seq.append(ligne);
    }
    preums = false;
  }

  m.insert(map<string, string>::value_type(sid, seq));

  if(verbose==true){
    fprintf(stderr, "done.\n");
    fflush(stderr);
  }
  
}



void loadFastaUltraFast(map<string, string> & m, char *fname, bool verbose){

  if(verbose==true){
    fprintf(stderr, "Loading %s ... ", fname);
    fflush(stderr);
  }

  string ligne, sid, seq;

  ifstream fic(fname);
  seq = "";
  while( getline(fic, ligne) ){
    sid = ligne.substr(1, SIZE_MAX_ID);
    getline(fic, seq);
    m.insert(map<string, string>::value_type(sid, seq));
  }

  if(verbose==true){
    fprintf(stderr, "done.\n");
    fflush(stderr);
  }
  
}




string my_reverse_cpt(string s, char *tab_cpt){

  string res = s;
  int lg = s.size();
  int i;

  for(i=0 ; i<lg ; i++){
    res[i] = tab_cpt[s[(lg-i)-1]];
  }

  return res;
}


void init_tab_cpt(char *tab_cpt){

  for(int i=0 ; i<'z' ; i++){
    tab_cpt[i] = 'N';
  }

  tab_cpt['a'] = 't';
  tab_cpt['t'] = 'a';
  tab_cpt['c'] = 'g';
  tab_cpt['g'] = 'c';

  tab_cpt['A'] = 'T';
  tab_cpt['T'] = 'A';
  tab_cpt['C'] = 'G';
  tab_cpt['G'] = 'C';

}


/*
loadFastaQ : loads a file in fastaq format (fasta + quality)
 */
int loadFastaQ(map<string, fastqRead> & m, char *fname, bool verbose){

  if(verbose==true){
    fprintf(stderr, "Loading %s ... ", fname);
    fflush(stderr);
  }
  
  string sid, seq, lsep, sq;
  string fullID, id, sep;
  long nbe = 0;
  
  ifstream fic(fname);
  while( getline(fic, sid) ){

    int posEnd = sid.find(' ');
    if( posEnd < 1 ){
      posEnd = sid.size();
    }

    fullID = sid.substr(1, posEnd-1);

    int posSep = fullID.find_first_of(' ');
    if(posSep==string::npos){ 
      posSep = fullID.find_first_of('\t');
    }
    if(posSep==string::npos){
      posSep = fullID.size();
    }
    id = fullID.substr(0, posSep);

    //fprintf(stderr, "fullID:'%s'\nID:'%s'\n", fullID.c_str(), id.c_str());

    getline(fic, seq);
    getline(fic, lsep);
    sep = lsep.substr(1, SIZE_MAX_ID);
    getline(fic, sq);
    
    //if(id!=sep){
      //fprintf(stderr, "ERROR in loadFastaq : %s != %s\nFasta read aborted.\n", id.c_str(), sep.c_str());
      //return 1;
      //}

    if(seq.size() != sq.size()){
      fprintf(stderr, "ERROR in LoadFastaQ : %s has size different from %s\n", seq.c_str(), sq.c_str());
      return 1;
    }

    fastqRead fr(id, fullID, seq, sq);

    m.insert(map<string, fastqRead>::value_type(id, fr));
    
    nbe++;
    if(verbose==true && nbe%1000000==0){
      fprintf(stderr, "%ld\n", nbe);
    }
  }

  if(verbose==true){
    fprintf(stderr, "done.\n");
    fflush(stderr);
  }

  return 0;
}


string my_reverse(string s){

  string res = s;
  int i, lg;

  lg = s.size();
  for(i=0 ; i<lg ; i++){
    res[i] = s[(lg-i)-1];
  }

  return res;
}


void myPrintFastaLine(string s, int NB_MAX, FILE *nf){

  int i, lg;
  lg = s.size();
  for(i=0 ; i<lg ; i+=NB_MAX){
    string stmp = s.substr(i, NB_MAX);
    fprintf(nf, "%s\n", stmp.c_str());
  }

}


void myPrintFastaMap(map<string,string> &mg, FILE *nf, int NB_MAX){

  map<string,string>::iterator it;

  for(it=mg.begin() ; it!=mg.end() ; it++){
    string chrID = it->first;
    string seq = it->second;
    fprintf(nf, ">%s\n", chrID.c_str());
    myPrintFastaLine(seq, NB_MAX, nf);
  }

}



fastqRead::fastqRead(void){
  ;
}


fastqRead::fastqRead(string &id, string &fullID, string &seq, string &squal){

  this->make(id, fullID, seq, squal);

}

void fastqRead::make(string &id, string &fullID, string &seq, string &squal){
  //istringstream iss(id_line);
  //string id;
  //getline(iss, id, ' ');

  this->fullID_ = fullID;

  this->id_ = id;
  this->seq_ = seq;
  this->squal_ = squal;

}



void fastqRead::reverse(char *tab_cpt){
  this->seq_ = my_reverse_cpt(this->seq_, tab_cpt);
  this->squal_ = my_reverse(this->squal_);
}

void fastqRead::fullPrint(FILE *nf){

  fprintf(nf, 
	  "@%s\n%s\n+\n%s\n",
	  this->fullID_.c_str(),
	  this->seq_.c_str(),
	  this->squal_.c_str()
	  );      
  
}


int fastqRead::trimLowQuality(char minQ){

  int lg = this->seq_.size();
  int posEnd = lg - 1;

  while( posEnd>=0 && this->squal_[posEnd] <= minQ ){
    posEnd--;
  }

  posEnd++;

  if( posEnd < lg ){
    this->squal_.erase(posEnd);
    this->seq_.erase(posEnd);
  }

  return (lg-posEnd);
}


bool fastqRead::readNext(ifstream &fic){

  string line;

  if( getline(fic, line) ){
    this->fullID_ = line;
    this->id_ = utils_fasta_extractID(line);
  } else {
    return false;
  }

  if( !getline(fic, this->seq_) ){ return false; }
  if( !getline(fic, line) ){ return false; }
  if( !getline(fic, this->squal_) ){ return false; }

  return true;
}


//-----------------------------------------------------------------------------------


double utils_fasta_getGC(string seq){

  long nbGC = 0;
  long nbAT = 0;

  long lg = seq.size();
  for(long ll=0 ; ll<lg ; ll++){
    char base = toupper(seq[ll]);
    if(base=='G' || base=='C'){
      nbGC++;
    } else {
      if(base=='A' || base=='T'){
	nbAT++;
      }
    }
  }

  long nb = nbGC + nbAT;

  double gc = (double)nbGC/(double)nb;

  return gc;

}


string utils_fasta_extractGeneID(string s, char delim){

  string stmp, geneID;
  
  istringstream iss(s);
  getline(iss, stmp, delim);
  int pp = stmp.find_first_of('.');
  if(pp==string::npos){ pp = stmp.size(); }
  geneID = stmp.substr(0, pp);
  
  return geneID;
 
}


string utils_fasta_extractTranscriptID(string s, char delim){
  string stmp, tID;

  if( s.find(delim)==string::npos ){
    return s;
  }

  istringstream iss(s);
  getline(iss, stmp, delim);
  getline(iss, tID, delim);

  return tID;
}


string utils_fasta_extractID(string &s){

  int posEnd = s.find(' ');
  if( posEnd==string::npos ){
    posEnd = s.find('\t');
  }

  if(posEnd==string::npos){
    posEnd = s.size();
  }

  string id = s.substr(1, posEnd-1);

  return id;
}


void loadGeneticCode(char *fName, map<string, string> &m){

  string line, codon, aa;
  ifstream fic(fName);
  while( getline(fic, line) ){
    istringstream iss(line);
    getline(iss, codon, '\t');
    getline(iss, aa, '\t');
    m[codon] = aa;
  }

}

#endif


