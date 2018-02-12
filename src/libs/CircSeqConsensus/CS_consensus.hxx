#ifndef _CS_consensus_hxx_
#define _CS_consensus_hxx_


/*
class CS_consensus stores the informaton on consensus sequences derived from repeats created by circle-seq sequencing.

Briefly, the consensus sequence is stored as a string and is derived from the RNAseq read sequence.
Consensus sequence can be derived from both paired-end and single-end reads.

For each position in the consensus sequence, I store all the corresponding base calls and quality scores (in the 'mc_' map). This information is 
used to derive the probability that a given base call in the consensus sequence is erroneous.

 */


#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>

#include "samRead/samRead.hxx"
#include "GFF/GFFtranscripts.hxx"


class CS_consensus{
public:
  CS_consensus(void);
  CS_consensus(const std::string &Read, const std::string &Qual, const int minRepeatSize, const int maxMM, const double minID);
  CS_consensus(const std::string &leftRead, const std::string &leftQual, const std::string &rightRead, const std::string &rightQual, const int minRepeatSize, const int maxMM, const double minID);
  CS_consensus(const std::string &read, const std::string &qual, const std::string csSeq, const double minID);
  CS_consensus(const std::string &leftRead, const std::string &leftQual, const std::string &rightRead, const std::string &rightQual, const std::string csSeq, const double minID);

  ~CS_consensus(void);
  void clean(void);

  bool make(const std::string &Read, const std::string &Qual, const int minRepeatSize, const int maxMM, const double minID);
  bool make(const std::string &leftRead, const std::string &leftQual, const std::string &rightRead, const std::string &rightQual, const int minRepeatSize, const int maxMM, const double minID);
  bool make(const std::string &read, const std::string &squal, const std::string &csSeq, const double minID);
  bool make(const std::string &leftRead, const std::string &leftQual, const std::string &rightRead, const std::string &rightQual, const std::string &csSeq, const double minID);
  bool make(const std::string &readSeq, const std::string &readQual, const std::string &mateSeq, const std::string &mateQual, const int csLg, const int breakPos, const int posCsInMate);
  bool make(const std::string &readSeq, const std::string &readQual, const int csLg, const int breakPos);

  void updateMap(const std::string &seq, const std::string &qual, const int posFirstRepeat, const int repeatSize);
  void updateConsensus(void);

  void prettyPrint(FILE *nf, const bool printQual);

  std::string getConsensus(void){ return this->consensus_; };
  std::string getCsQuality(const int minQual, const int maxQual);
  int getCsProb(std::vector<int> &vQual);
  int getCsProb(std::vector<int> &vQual, std::vector<int> &vDepth, std::vector<int> &vDivQual);

  int getPosStartInRead(void);
  int getPosStartInMate(void);

  void reverseCpt(const char *tab_cpt);

  int findFirstMismatch(std::map<int,int> &mc, std::string &chrSeq);
  int moveToRight(std::map<int,int> &mTmp, std::map<int,int> &mc,  int nbToMove, GFFtranscript & gt);
  int moveToLeft(std::map<int,int> &mTmp, std::map<int,int> &mc, int nbToMove, GFFtranscript & gt);

  int refine(int *offset, samRead &sr, vector<GFFtranscript> &vt, map<string,string> &mg);
  //int refine(std::string &upSeq, std::string &middleSeq, std::string &downSeq, std::string &scs);
  //int refine(std::map<int,int> &mc, int posMMinCs, std::string &chrSeq, std::map<int,int> &refinedMap, std::string &newCs);
  //int refine(std::map<int,int> &mc, std::string &chrSeq, int *offset);
  //int refineWithIndels(std::map<int,int> &mc, std::string &chrSeq, int *offset, samRead & sr);
  //int refine(std::map<int,int> &mc, std::string &chrSeq, int *offset, samRead &sr, GFFtranscript &gt);

  int getSize(void);

  double getPrctId(void);
  void setLimitHomo(void);
  void setLimitHomo(int left, int right);
  void getHomoSizes(int *left, int *right);

private:
  std::pair<int,int> limitHomo_;
  std::pair<int,int> posStartFirst_;
  std::string consensus_;
  std::map<int, std::pair<std::string, std::string> > mc_; /**< a map containing the list of basecalls and qualities at each position in the consensus */
};


#endif
