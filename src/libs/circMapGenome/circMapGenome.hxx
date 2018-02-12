#ifndef _CIRCMAPGENOME_HXX_
#define _CIRCMAPGENOME_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>

#include "UTILS_FASTA/utils_fasta.hxx"
#include "GFF/mapOfStrand.hxx"

#include "circMapCov-simple/circMapCovSimple.hxx"
#include "circMapIndels/circMapIndels.hxx"
#include "circMapObs/circMapObs.hxx"

using namespace std;

class circMapGenomePos{
public:
  circMapGenomePos(void);
  circMapGenomePos(simplePosCov &spc);
  circMapGenomePos(C_indelPos &ip);
  circMapGenomePos(long nbObs);

  void init(void);
  void init(simplePosCov &spc);
  void init(C_indelPos &ip);
  void init(long nbObs);

  char initRefBase(map<string,string> &mg, string chrID, long gPos, mapOfStrand &mos, char *tab_cpt);
  void setRefBase(char refBase);

  void set(simplePosCov &spc);
  void set(C_indelPos &ip);
  void set(long nbObs);

  void add(C_indelPos &ip);
  void add(simplePosCov &spc);
  void add(long nbObs);

  long getNbObservations(void);

  //private:
  simplePosCov spc_;
  C_indelPos ip_;
  long nbObs_;
};


class circMapGenome{
public:
  circMapGenome(void);
  circMapGenome(char *fname);
  circMapGenome(map<string,string> &mg);

  void init(void);

  // Map of coverage functions
  void loadMapCoverage(char *fname, bool add);
  void printMapCoverage(FILE *nf, const char *tabBases);
  int updateMapCoverage(const char *fileSam, const int *tabCor, const int minScore);
  int processSamRead_forCoverage(const samRead &sr, const int* tabCor, const int minScore);
  int processSamRead_forCoverage(const samRead &sr, map<string,map<long,circMapGenomePos> >::iterator itc, const int* tabCor, const int minScore);
  void getRefAndAlterCounts_cov(const std::string &chrID, const int gPos, int *nbRef, int *nbAlter, const char *tabBases);
  void getRefAndAlterCounts_cov(map<string, map<long,circMapGenomePos> >::iterator itc, const int gPos, int *nbRef, int *nbAlter, const char *tabBases);
  bool isSuspect_cov(const std::string &chrID, const int gPos, const int minCov, const double maxFrac, const char *tabBases, int *nbRef, int *nbTot);
  bool isSuspect_cov(map<string, map<long,circMapGenomePos> >::iterator it, const int gPos, const int minCov, const double maxFrac, const char *tabBases, int *nbRef, int *nbTot);

  bool updateAllRefBase(map<string,string> &mg, mapOfStrand &mos, char *tab_cpt);
  bool updateAllRefBase(map<string,string> &mg);

  // Map of indels functions
  void loadMapIndels(char *fname, bool add=false);
  void printMapIndels(FILE *nf);
  int updateMapIndels(char *fname);
  int updateIndels(string &chrID, map<int,pair<int, int> > &mIns, map<int, pair<int,int> > &mDel);
  int updateIndels(map<string, map<long, circMapGenomePos> >::iterator itc, map<int,pair<int, int> > &mIns, map<int, pair<int,int> > &mDel);
  int updateMapIndels(samRead &sr);
  int updateMapIndels(samRead &sr, map<string, map<long,circMapGenomePos> >::iterator itc);
  int updateCovIndels(string &chrID, map<int,int> &mc);
  int updateCovIndels(map<string, map<long,circMapGenomePos> >::iterator itc, map<int,int> &mc);
  bool isSuspect_indels(std::string chrID, const long gPos, const long minCoverage, const double maxFracAlter);
  bool isSuspect_indels(map<string, map<long,circMapGenomePos> >::iterator itm, const long gPos, const long minCoverage, const double maxFracAlter);

  // Map of observations functions
  void loadObservations(char *fname, int nbLinesOfHeader=0, bool erasePreviousMap=0);
  void addObservations(char *fname, int nbLinesOfHeader);
  void printObservations(FILE *nf, bool printZeros=false);
  void printObservationsWithBase(FILE *nf, map<string,string> &mg, mapOfStrand &mos, bool printStrand);
  int addObservation(std::string chrID, long gPos);
  int addObservation(map<string, map<long,circMapGenomePos> >::iterator itc, long gPos);
  int addObservation(std::string chrID, vector<long> & vpos);
  int addObservation(map< string,map<long,circMapGenomePos> >::iterator itc, vector<long> & vpos);

  //private:
  map<string, map<long, circMapGenomePos> > mg_;
};




#endif
