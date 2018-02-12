
#ifndef _UTILS_FASTA_HXX_
#define _UTILS_FASTA_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>

#define SIZE_MAX_ID 200

using namespace std;


class fastqRead{

public:
  fastqRead(void);
  fastqRead(std::string &id, std::string &fullID, std::string &seq, std::string &squal);

  void make(std::string &id, std::string &fullID, std::string &seq, std::string &squal);

  void fullPrint(FILE *nf);
  void reverse(char *tab_cpt);
  bool readNext(ifstream &fic);
  int trimLowQuality(char minQ);

  std::string fullID_;
  std::string id_;
  std::string seq_;
  std::string squal_;

};



void loadFasta(map<string, string> & m, char *fname, bool extractID=true, bool verbose=false, bool debug=false);
void loadFasta(map<string, string> & m, char *fname, string chrID, bool extractID=true, bool verbose=false, bool debug=false);
void loadFastaUltraFast(map<string, string> & m, char *fname, bool verbose);
int loadFastaQ(map<string, fastqRead> & m, char *fname, bool verbose);

string my_reverse_cpt(string s, char *tab_cpt);
string my_reverse(string s);

void init_tab_cpt(char *tab_cpt);

double utils_fasta_getGC(string seq);

void myPrintFastaLine(string s, int NB_MAX, FILE *nf);
void myPrintFastaMap(map<string,string> &mg, FILE *nf, int NB_MAX=60);

string utils_fasta_extractGeneID(string s, char delim);
string utils_fasta_extractTranscriptID(string s, char delim);
string utils_fasta_extractID(string &s);


void loadGeneticCode(char *fName, map<string, string> &m);

#endif
