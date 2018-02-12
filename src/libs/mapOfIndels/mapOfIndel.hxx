#ifndef _mapOfIndel_hxx_
#define _mapOfIndel_hxx_


#include <vector>
#include <map>
#include <iostream>

#include "indelCoverage.hxx"
#include "remappedRead/remappedRead.hxx"

using namespace std;


class C_mapOfIndel{

public:
  C_mapOfIndel(void);
  C_mapOfIndel(map<string, string> &mg, bool verbose);

  void init(map<string, string> &mg, bool verbose);


  //void updateFromRemapped(char *fname, map<string, string> &mBarcodes, bool verbose);
  //void updateFromRemapped(map<string, vector< pair<remappedRead, remappedRead> > > & mr, map<string, string> &mBarcodes, bool verbose);
  void updateFromFamily( vector< pair<remappedRead, remappedRead> > &v);
  void updateFromSamRead( pair<samRead, samRead> &pp );
  int updateFromFamily(
		       vector< vector< pair<remappedRead, remappedRead> > > &vFam,
  		       FILE *nfCandidates,
  		       long *NB_INS,
  		       long *NB_DEL,
  		       long *NB_POS_COVERED_FOR_INDELS,
  		       map<string, string> &mg,
  		       char *tab_cpt,
  		       char *tabBases
  		       );



  void updateFromPreviousMap(char *fname);

  void display(FILE *nf);

  map<string, vector<C_indelCoverage> > indelCoverage_;
};





//void initMapConsensusIndel(map<long, char> &m, remappedRead &r);
//void updateMapConsensusIndel(map<long, char> &m, remappedRead &r);
void getMapConsensusIndel(map<long, char> &m, vector< pair<remappedRead, remappedRead> > &vFamily);




int updateIndelErrors(
		       C_mapOfIndel &moi,
		       int minCov,
		       double maxProp,
		       vector < vector< pair<remappedRead, remappedRead> > > &vFam,
		       FILE *nfCandidates,
		       long *NB_INS,
		       long *NB_DEL,
		       long *NB_POS_COVERED_FOR_INDELS,
		       map<string, vector<int> > &mhp, // map of homopolymer size for every position in every transcript (key=transcriptID)
		       map<int, pair<long,long> > &mObsHp, // map of observations for every size of homopolymer run
		       map<string, GFFtranscript> &mGFF, // map of transcripts (needed to convert genomic to transcript coordinates)
		       map<string, vector<long> > &mObs,
		       map<string, string> &mg,
		       char *tab_cpt,
		       char *tabBases
		      );


void displayIndelCandidate(FILE *nf, vector< vector< pair<remappedRead, remappedRead> > > &vFam, char stat, int indelSize, string indel, string chrID, long gPos, string tID, long tPos, int runSize, string strand, long cov, long ins, long del, double prop);

void displayIndelCandidate(FILE *nf, vector< vector< pair<remappedRead, remappedRead> > >&vFam, char stat, int indelSize, string chrID, long gPos, string tID, long tPos, int runSize, string strand, long readCoord, long mateCoord);

int getListOfIndels(
		     vector< pair<remappedRead, remappedRead> > &vr, 
		     vector< pair<long, long> > &vIns,
		     vector< pair<long, long> > &vDel,
		     bool onlyOverlap
		    );

void getSharedIndels(vector< vector< pair<remappedRead, remappedRead> > > &vFam, 
		     int NB_ROUND,
		     vector< pair<long, long> > &vIns, 
		     vector< pair<long, long> > &vDel,
		     pair<remappedRead, remappedRead> &prr
		     );

void initIndelsTmp(pair<remappedRead, remappedRead> &prr, vector< pair<long, long> > &vIns, vector< pair<long, long> > &vDel);
void updateIndelsTmp(pair<remappedRead, remappedRead> &prr, vector< pair<long, long> > &vIns, vector< pair<long, long> > &vDel);

void displayMSC(FILE *nf, map<long, char> &msc, map<string, string> &mg, string &chrID, bool cpt, char *tab_cpt);

#endif
