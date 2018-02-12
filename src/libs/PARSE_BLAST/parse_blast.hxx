#ifndef _PARSE_BLAST_HXX_
#define _PARSE_BLAST_HXX_

#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>
#include <algorithm>


class blast_m8{
public:
  blast_m8();
  blast_m8(std::string &s);

  void make(std::string &s);
  void print(FILE *nf);

  std::string queryId_;
  std::string subjectId_;
  double prctIdent_;
  int alnLength_;
  int mismatchCount_;
  int gapOpenCount_;
  int queryStart_;
  int queryEnd_;
  int subjectStart_;
  int subjectEnd_;
  std::string strEval_;
  double eVal_;
  double bitScore_;
};


#endif
