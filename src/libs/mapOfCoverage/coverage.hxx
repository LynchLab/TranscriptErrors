#ifndef _coverage_HXX_
#define _coverage_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>

//#include "boost/archive/binary_iarchive.hpp"
//#include "boost/archive/binary_oarchive.hpp"
//#include "boost/serialization/map.hpp"
//#include "boost/serialization/vector.hpp"


using namespace std;
//using namespace boost;
//using namespace serialization;

#define _STRAND_PLUS_ '+'
#define _STRAND_MINUS_ '-'
#define _STRAND_UNKNOWN_ '?'
#define _STRAND_BOTH_ 'b'


#define NB_BASES_MOC 4

class coverage{
public:
  coverage(void);
  coverage(long val);
  coverage(long *tabVal);
  //char strand;
  //char genomicBase;

  void add(long *tabVal);
  void add(char c, char *tabBases);
  void display(FILE* nf);

  void addFromGenomic(coverage &cv, char *tabBases, char strand, char *tabCpt);

  long adjustedCov_[NB_BASES_MOC];

  //template<typename Archive>
  //void serialize(Archive& ar, const unsigned version) {
  //ar & adjustedCov_; 
  //}

};



#endif
