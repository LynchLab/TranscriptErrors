#ifndef _CIRC_MAP_OBS_
#define _CIRC_MAP_OBS_

/*

  the circMapObs class is used to store the 'observations' in a circle-seq experiment.
  1 observation = 1 usable alignment between a consensus sequence and the reference genome (usable = that passes the criteria for consensus quality / heterozygocity / ...)

  For each position in the genome, I record the number of observations.

  The data is stored in a map of vector. Each entry in the map corresponds to a chromosome (scaffold).

  On disk, the map is printed in a 3 columns files with tab separation (chromosomeID position #observations)

 */


#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>

#include "UTILS_FASTA/utils_fasta.hxx"
#include "GFF/mapOfStrand.hxx"


class circMapObs{
public:
  circMapObs(void);
  circMapObs(std::map<std::string,std::string> &mg);
  circMapObs(std::vector<std::string> & vFiles);
  circMapObs(char *fname);

  void init(std::map<std::string,std::string> &mg);
  void init(std::vector<std::string> & vFiles);
  void init(char *fname);
  void add(char *fname);


  void addObs(std::string chrID, int pos);
  void addObs(std::string chrID, std::vector<long> vpos);
  void setObs(std::string chrID, int pos, long nb);

  long getObs(std::string chrID, int pos);

  long getSum(void);
  void print(FILE *nf);
  void printWithBase(FILE *nf, std::map<std::string,std::string> & mg, mapOfStrand &mos);

private:
  std::map<std::string, std::vector<long> > mObs_;

};


class circMapObsAndBase{
public:
  circMapObsAndBase(void);
  circMapObsAndBase(char *fname);

  void init(char *fname);
  int get(std::string chrID, long gPos, std::pair<long,char> &pp);

private:
  std::map<std::string , std::map<long, std::pair<long,char> > > m_;
};

#endif
