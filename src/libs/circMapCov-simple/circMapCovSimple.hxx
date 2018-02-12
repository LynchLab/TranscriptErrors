#ifndef _CIRC_MAP_COV_SIMPLE_HXX_
#define _CIRC_MAP_COV_SIMPLE_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>

#include "UTILS_FASTA/utils_fasta.hxx"
#include "samRead/samRead.hxx"


class simplePosCov{
public:
  simplePosCov(void);
  simplePosCov(const char base);

  int init(const char base);

  int add(const int iBase);
  int add(const char base, const int *tabCor);

  int getCoverage(const int iBase);
  int getCoverage(const char base, const int *tabCor);

  char getRefBase(void);
  void setRefBase(char refBase);

  void getRefAndAlterCounts(int *nbRef, int *nbAlter, const char* tabBases);

  bool isSuspect(const int minCov,const double maxFrac, const char *tabBases, int *nbRef, int *nbTot);

  void print(FILE *nf) const;

  int set(const int iBase, const int nb);
  int set(const char base, const int nb, const int*tabCor);

  void add(simplePosCov &spc);

private:
  char refBase_;
  int cov_[4];
};


class circMapCovSimple{
public:
  circMapCovSimple(void);
  circMapCovSimple(char *fileGenome);
  circMapCovSimple(const std::map<std::string,std::string> &mg);
  circMapCovSimple(char *fileGenome, char *fileSam, const int *tabCor, const int minScore=30);

  void load(char *fname);
  void load(char *fname, string chrOnly);

  int init(const std::map<std::string,std::string> &mg);

  int update(const char *fileSam, const int *tabCor, const int minScore=30);
  int processSamRead(const samRead &sr, const int* tabCor, const int minScore=30);

  void getRefAndAlterCounts(const std::string &chrID, const int gPos, int *nbRef, int *nbAlter, const char *tabBases);
  void getRefAndAlterCounts(std::map<std::string, std::map<long,simplePosCov> >::iterator it, const int gPos, int *nbRef, int *nbAlter, const char *tabBases);

  bool isSuspect(const std::string &chrID, const int gPos, const int minCov, const double maxFrac, const char *tabBases, int *nbRef, int *nbTot);
  bool isSuspect(std::map<std::string, std::map<long,simplePosCov> >::iterator it, const int gPos, const int minCov, const double maxFrac, const char *tabBases, int *nbRef, int *nbTot);
  bool isSuspect(const std::map<std::string , std::map<long,simplePosCov> >::const_iterator it, const int gPos, const int minCov, const double maxFrac, const char *tabBases, int *nbRef, int *nbTot);


  void print(FILE *nf, const char *tabBases) const;

private:
  std::map<std::string, std::map<long,simplePosCov> > mc_;
};


#endif
