#ifndef _mapOfCOverage_HXX_
#define _mapOfCOverage_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>


using namespace std;

#include <stdlib.h>

#include "coverage.hxx"
#include "remappedRead/remappedRead.hxx"
#include "remappedRead/remappedRead_utils.hxx"
#include "samRead/samRead.hxx"


class mapOfCoverage{
public:
  mapOfCoverage(void);

  void init(map<string, string> &mg, bool verbose);

  void updateFromRemapped(char *fname, map<string, string> &mg, char *tab_cpt, map<string, string> &mBarcodes, char *tabBases, bool verbose);
  void updateFromRemapped(map<string, vector< pair<remappedRead, remappedRead> > > & mr, map<string, string> &mg, char *tab_cpt, map<string, string> &mBarcodes, char *tabBases, bool verbose);
  void updateFromPreviousMap(char *fname);
  void updateFromRemapped(int nbRounds, map<string, vector< pair<remappedRead, remappedRead> > > * mr, map<string, string> &mg, char *tab_cpt, map<string, string> * mBarcodes, char *tabBases, bool verbose);
  void updateFromFamily( vector< pair<remappedRead, remappedRead> > &vf, char *tab_cpt, char *tabBases);

  void update(map<long, char> &mc, string &chrID, char *tab_cpt, bool cpt, char *tabBases);
  void updateFromNoGapSamRead(pair<samRead, samRead> &pp, char *tabBases, int SHIFT);


  void mergeTOgDNAmap(char *genomicMapFile, map<string, GFFtranscript> &mst, map<string, string> &mg, char *tabBases, char *tab_cpt);
  //void mergeTOgDNAmap(map<string, vector<coverage> > &mgd, map<string, vector<char> > &mStrand, char *tabBases, char *tab_cpt);
  void mergeTOgDNAmap(mapOfCoverage &mGenomic, map<string, vector<char> > &mStrand, char *tabBases, char *tab_cpt);


  void getCoverageAt(string chrID, long pos, coverage &c);

  void display(FILE *nf, char *tabBases, map<string, string> &mg);
  void display(FILE *nf, char *tabBases);
  void displayHeader(FILE *nf, char *tabBases);

  void getCovAndProp(string &chrID, long pos, bool cpt, char gc, char *tabBases, long *totalCov, double *prop);

  map<string, vector<coverage> > coverage_;


  //////////////////////////////////////////////////////////////////
  // BOOST serialization functions:
  //template<typename Archive>
  //void serialize(Archive& ar, const unsigned version) {
  //ar & coverage_;
  //}

  //void binarySave(char *fname);
  //void binaryLoad(char *fname);

};


void init_mStrand(map<string, vector<char> > &mStrand, map<string, vector<coverage> > &mCov);
void load_mStrand(map<string, vector<char> > &mStrand, map<string, GFFtranscript> &mst);

//void binaryLoad(char *fname, mapOfCoverage *m);

#endif
