#ifndef _GFFTRANSCRIPT_HXX_
#define _GFFTRANSCRIPT_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>

#include "GFFentry.hxx"


using namespace std;

class GFFtranscript{

public:
  GFFtranscript(void);
  GFFtranscript(vector<GFFentry> & vGFF, string sTAG, int shift, long posStart=-1, long posEnd=-1);
  GFFtranscript(vector<GFFentry> & vGFF, string sTAG, int shift, map<string, string> & mg, char *tab_cpt, long posStart=-1, long posEnd=-1);
  
  void init(vector<GFFentry> & vGFF, string sTAG, int shift, map<string, string> & mg, char *tab_cpt, long posStart=-1, long posEnd=-1);

  void getAllEEJcoords(vector<long> &vCoords);
  void getAllIntronsCoordsInUnSplicedTranscripts(vector< pair<long, long> > &vCoords);

  long getSplicedSize(void);

  long genomicToSplicedTranscript(long pos);
  long genomicToUnSplicedTranscript(long pos); // KNOWN BOG: does not work for position within intron

  long splicedToGenomic_deprecated(long pos);
  long splicedToGenomic(long pos);
  int allSplicedToGenomic(map<int,int> &m);

  string getUnSplicedSeq(map<string, string> &mg, char *tab_cpt);
  string getSplicedSeq(void);

  string getIntron(int intronNumber, map<string, string> &mg, char *tab_cpt);
  void getAllIntrons(vector<string> & vIntrons, map<string, string> &mg, char *tab_cpt);

  string getExon(int exonNumber, map<string, string> &mg, char *tab_cpt);
  void getAllExons(vector<string> &vExons, map<string, string> &mg, char *tab_cpt);

  string getIntron(int intronNumber); // TO DO !!!
  void getAllIntrons(vector<string> & vIntrons); // TO DO !!!

  void getSpliceUnSpliceAndFlanking(int intronNumber, int flankingSize, map<string, string> &mg, char *tab_cpt, string &spliced, string &unSpliced);
  void getFlankingFromIntron(int intronNumber, int flankingSize, map<string, string> &mg, char *tab_cpt, string &flank5p, string &flank3p);

  string getIntronWithFlanking(int intronNumber, int flankingSize, map<string, string> &mg, char *tab_cpt);
  bool intronIsNeverExonic(int intronNumber, map<string, vector<bool> > &mpi, bool exonic);

  void markIntronsGenomicScale(vector<bool> &v, bool mark);
  void markExonsGenomicScale(vector<bool> &v, bool mark);

  void getGenomicMinMaxPos(long *start, long *end);
  string getUpStreamGenomicSeq(int nbBases, string &chromosomeSeq, char *tab_cpt, int *nbDone);
  string getDownStreamGenomicSeq(int nbBases, string &chromosomeSeq, char *tab_cpt, int *nbDone);


  int getNbExons(void);
  int getNbIntrons(void);
  bool isOnPlusStrand(void);

  void printCoords(FILE *nf);

  
  string id_;
  string source_;
  string chr_;
  string strand_;
  string seq_;
  vector< pair<long, long> > vexonsCoords_;

  long posMin_, posMax_;  
  long posStart_, posStop_;
};



bool my_GFFentry_compare(GFFentry e1, GFFentry e2);


#endif
