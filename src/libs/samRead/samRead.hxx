#ifndef _SAM_READ_HXX_
#define _SAM_READ_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>

#include "GFF/GFFtranscripts.hxx"

void myLocal_samRead_printMap(FILE *nf, map<int,int> & m);

class samRead{

public:
  samRead(void);
  samRead(std::string line);
  ~samRead(void);

  int init(std::string line);
  bool initWithNext(std::ifstream &f);
  int readNext(std::ifstream &fic);
  int readNext(std::ifstream &fic, std::string &line);


  bool getGenomicCoords(std::map<int,int> &m) const; // potential bog, function retired
  std::string getGenomicSeq(std::map<int,int> &mc, std::map<std::string , std::string> &mg) const;

  bool isPaired(void);
  bool readMappedInProperPair(void);
  bool readUnmapped(void);
  bool mateUnmapped(void);
  bool readReverseStrand(void);
  bool mateReverseStrand(void);
  bool firstInPair(void);
  bool secondInPair(void);
  bool notPrimaryAlignment(void);
  bool readFailsQualityChecks(void);
  bool readIsPCRduplicate(void);
  bool supplementaryAlignment(void);

  bool isOnIntron(void);

  std::string getQueryName(void) const;
  std::string getRefName(void) const;
  std::string getCigar(void) const;
  std::string getSeq(void) const;
  int getSeqLg(void) const;
  int getFlag(void) const;
  int getNbMM(void);
  int getEditDistanceFromFlag(void);
  int getEditDistance(void);
  int getEditDistance(int *nbMM, int *nbIndels, int *nbIns, int *nbDel, int *totalSizeIndels, int *totalSizeIns, int *totalSizeDel, int *nbN, int *nbOther);
  int getPos(void);
  void setPos(int pos);

  bool getSplitCigar(std::vector< std::pair<int, char> > & v) const;
  int getCigarDepth(void);

  //int getIndelPosAndSize(int *pos, int *size);
  //int getIndelsCoords(std::map<int,int> &mIns, std::map<int,int> &mDel) const;
  int getEditDistance(map<int,int> & m, string & chrSeq);
  int getNbIndels(void);
  int getIndels(std::map<int, std::pair<int,int> > &mIns, std::map<int, std::pair<int,int> > &mDel);
  int getIndels(std::map<int, std::pair<int,int> > &mIns, std::map<int, std::pair<int,int> > &mDel, std::map<int,int> &mc);

  int refine(std::map<int,int> &m, string & chrSeq, GFFtranscript &gt, int *offset);
  int findFirstMismatch(std::map<int,int> &m, std::string &chrSeq);
  int findFirstEdit(std::map<int,int> &m, std::string &chrSeq);

  int getMapOfCoords(std::map<int,int> &m) const;
  int getMapOfCoordsV2(std::map<int,int> &m);

  int moveToLeft(std::map<int,int> &m, std::map<int,int> &newM, int * nbToMove, GFFtranscript &gt, std::string &chrSeq, int indelType, int indelPos, int indelSize);
  int moveToRight(std::map<int,int> &m, std::map<int,int> &newM, int * nbToMove, GFFtranscript &gt, std::string &chrSeq, int indelType, int indelPos, int indelSize);


  int tryToFillDeletionLeft(std::map<int, int> &newM, int indelPos, int indelSize, GFFtranscript &gt, string &chrSeq);

  int printAln(FILE *nf, std::string &chrSeq);
  int printAln(FILE *nf, std::map<int,int> &m, std::string &chrSeq);

  std::map<std::string,std::string> mOpt_;

private:
  std::string qName_;
  int flag_;
  std::string rName_;
  int pos_;
  int mapq_;
  std::string cigar_;
  std::string rNext_;
  int pNext_;
  int tLen_;
  std::string seq_;
  std::string qual_;

};


#endif
