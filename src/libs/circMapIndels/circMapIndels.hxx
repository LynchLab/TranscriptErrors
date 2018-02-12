#ifndef _CIRC_MAP_INDELS_HXX_
#define _CIRC_MAP_INDELS_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>

#include "UTILS_FASTA/utils_fasta.hxx"
#include "samRead/samRead.hxx"


class C_indelPos{
public:
  C_indelPos(void);
  C_indelPos(const long cov, const long nbIns, const long nbDel);

  void get(long *cov, long *nbIns, long *nbDel);
  void set(const long cov, const long nbIns, const long nbDel);

  void add(C_indelPos &ip);
  void add(long cov, long nbIns, long nbDel);

  void addCov(void);
  void addCov(const long nb);
  void addIns(void);
  void addIns(const long nb);
  void addDel(void);
  void addDel(const long nb);

  bool isSuspect(const long minCovergae, const double maxFracAlter);

  void print(FILE *nf);

private:
  long cov_;
  long nbIns_;
  long nbDel_;
};


/*
class C_circMapIndels{
public:
  C_circMapIndels(void);
  C_circMapIndels(const char *genomeFile);
  C_circMapIndels(const map<string,string> &mg);


  int init(const char *genomeFile);
  int init(const std::map<std::string,std::string> &mg);

  int update(const char *samFile);


  int updateIndels(std::string &chrID, std::map<int, std::pair<int, int> > &mIns, std::map<int, std::pair<int,int> > &mDel);
  int updateIndels(std::map<std::string , std::map<long,C_indelPos> >::iterator itm , std::map<int, std::pair<int, int> > &mIns, std::map<int, std::pair<int,int> > &mDel);


  int updateCov(std::string &chrID, std::map<int,int> &mc);
  int updateCov(std::map<std::string, std::map<long,C_indelPos> >::iterator itm, std::map<int,int> &mc);


  void print(FILE *nf) const;
  void load(char *fname);
  void load(char *fname, string chrID);

  int addCov(const long gPos);
  int set(const std::string chrID, const long gPos, const long cov, const long nbIns, const long nbDel);
  int set(std::map<std::string, std::map<long,C_indelPos> >::iterator itm, const long gPos, const long cov, const long nbIns, const long nbDel);


  int get(const std::string chrID, const long gPos, long * cov, long * nbIns, long * nbDel);
  int get(std::map<std::string, std::map<long,C_indelPos> >::iterator itm, const long gPos, long * cov, long * nbIns, long * nbDel);

  bool isSuspect(std::string chrID, const long gPos, const long minCovergae, const double maxFracAlter);
  bool isSuspect(std::map<std::string, std::map<long,C_indelPos> >::iterator itm, const long gPos, const long minCovergae, const double maxFracAlter);

private:
  std::map<std::string, std::map<long,C_indelPos> > m_;
};
*/

#endif
