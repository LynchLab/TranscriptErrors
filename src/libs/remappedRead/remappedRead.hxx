#ifndef _REMAPPED_READ_HXX_
#define _REMAPPED_READ_HXX_

#include <stdlib.h>

#include <map>
#include <vector>
#include <algorithm>

#include <iostream>
#include <vector>


class remappedRead{

public:
  remappedRead(void);
  remappedRead(std::string, bool includeSequence=true, bool computeSeqWithIndels=false, bool includeQuality=false);

  void quickDisplay(FILE *nf);
  void fullPrint(FILE *nf, bool printQuality=false);
  bool coversPosition(long pos);
  long genomeToReadCoords(long pos);
  bool asInsertion(std::pair<long, long> &pi);
  bool asDeletion(std::pair<long, long> &pi);

  std::string getSeqWithIndels(bool qual=false);

  std::string qName_; // = read ID
  std::string strand_;
  std::string rName_; // Typically the chromosome/scaffold name
  std::string transcriptID_; // Just to keep track of the gene/transcript onto which the read maps
  //std::vector< std::pair<long,long> > qmaps_; // Vector of bloc coordinates on the read
  //std::vector< std::pair<long,long> > tmaps_; // Vector of corresponding bloc coordinates on the target
  std::vector< std::pair<long,long> > vmaps_;

  std::string seq_;
  std::string seqWithIndels_;

  std::string squal_;

  std::vector<std::pair<long,long> > vInsertions_;
  std::vector<std::pair<long,long> > vDeletions_;
};


#endif
