#ifndef _mapOfStrand_HXX_
#define _mapOfStrand_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>


#include "UTILS_FASTA/utils_fasta.hxx"
#include "GFF/GFFtranscripts_utils.hxx"


using namespace std;



class mapOfStrand{
public:
  mapOfStrand(void);
  mapOfStrand(map<string, string> &mg, map<string, GFFtranscript> &mt);

  void init(map<string, string> &mg);
  char getStrand(string chrID, long pos);

  map<string, vector<char> > m_;
};



#endif
