#ifndef _TRerrors_HXX_
#define _TRerrors_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>


#include "GFF/GFFtranscripts.hxx"
#include "remappedRead/remappedRead.hxx"
#include "mapOfCoverage/mapOfCoverage.hxx"



using namespace std;

class mapOfCoverage;




//void getConsensus(map<long, pair<char, char> > &mi, map<long, char> &mc, bool onlyOverlap);
//void getConsensus(remappedRead &read, remappedRead &mate, map<long, char> &mc, char *tab_cpt, bool onlyOverlap);
//void getConsensus(vector< map<long, char> > &vmc, map<long, char> &mc);
//void getConsensus(vector< pair<remappedRead, remappedRead> > &vpr, map<long, char> &mc, char *tab_cpt, bool onlyOverlap);


void getConsensus(vector< vector< pair<remappedRead, remappedRead> > > &vFam, map<long, char> &mhc, map<long, char> &msc, char *tabBases, char *tab_cpt, bool onlyOverlap);

void getConsensus(vector< pair<remappedRead, remappedRead> > &vr, map<long, char> &mhc, map<long, char> &msc, char *tabBases, char *tab_cpt, bool onlyOverlap);
void getExtendedConsensus(vector< pair<remappedRead, remappedRead> > &vr, map<long, string> &mc, char *tabBases, char *tab_cpt, bool onlyOverlap);
void getCollapsedConsensus(vector< pair<remappedRead, remappedRead> > &vr, map<long, string> &mc, char *tabBases, char *tab_cpt);

//void getSoftConsensus(vector< pair<remappedRead, remappedRead> > &vr, map<long, char> &msc, char *tabBases, char *tab_cpt, bool onlyOverlap);

void getHardAndSoftConsensus(string &s, char *hc, char *sc, char *tabBases);

//void prepareConsensus(map<long, char> &mc1, map<long, char> &mc);
//void addConsensus(map<long, char> &mcs, map<long, char> &mc);


int getIntersection(remappedRead &read, remappedRead &mate, map<long, pair<char, char> > &mi, char *tab_cpt);


long updateTranscriptionErrors(
			       long minCov,
			       double minProp,
			       bool onlyOverlap, 
			       map<string, vector< vector< pair<remappedRead, remappedRead> > > >::iterator itf, 
			       int NB_ROUNDS, 
			       map<string, GFFtranscript> &mt,
			       map<string, string> &mg, 
			       char *tab_cpt, 
			       mapOfCoverage &mcov,
			       char *tabBase,
			       FILE *nf_candidates, 
			       long *NB_TR_ERRORS, 
			       long *NB_POS_COVERED_FOR_TR_ERRORS,
			       map<string, vector<long> > &mObs
			       );


long updateTranscriptionErrors(
			       bool onlyOverlap,
			       long minCov,
			       double minProp,
			       map<long, char> &mhc,
			       map<long, char> &msc,
			       map<string, vector< vector< pair<remappedRead, remappedRead> > > >::iterator itf, 
			       map<string, GFFtranscript> &mt,
			       map<string, string> &mg,
			       char *tab_cpt,
			       mapOfCoverage &mcov,
			       char *tabBases,
			       FILE *nf_candidates,
			       long *NB_TR_ERRORS,
			       long *NB_POS_COVERED_FOR_TR_ERRORS,
			       map<string, vector<long> > &mObs
			       );



void debogConsensus(vector< pair<remappedRead, remappedRead> > &vFamily, map<long, char> &mhc, map<long, char> &msc);

/*
void updateTranscriptionErrors(
			       bool onlyOverlap,
			       map<string, vector< pair<remappedRead, remappedRead> > > * mNbBarcodes, 
			       string barcode,
			       int NB_ROUNDS, 
			       bool *roundHasPos, 
			       map<string, string> &mg,
			       char *tab_cpt,
			       FILE *nf_candidates, 
			       long *NB_TR_ERRORS, 
			       long *NB_POS_COVERED_FOR_TR_ERRORS
			       );
*/


/*
void updateMapFamilies(
		       int iRound,
			map<string, vector< vector< pair< remappedRead, remappedRead> > > > & mf,
			map<string, vector< pair<remappedRead, remappedRead> > > & mBarcodesAtPosition,
			vector<string> & vBarcodesSeen,
			int nbRounds
		       );
*/

/*
void displayCandidate(
		      FILE *nf,
		      map<string, vector< pair<remappedRead, remappedRead> > > *mNbBarcodes, 
		      string barcode, 
		      bool *roundHasPos,
		      int nbRounds,
		      string chrID,
		      long pos, 
		      char gc, 
		      char c);

*/

void displayCandidate(
		      FILE *nf_candidates, 
		      char gc, 
		      char c, 
		      long cov,
		      double prop,
		      string chrID, 
		      long gPos, 
		      string tID, 
		      long tPos, 
		      map<string, vector< vector< pair<remappedRead, remappedRead> > > >::iterator  itf,
		      long posInMap,
		      long mapSize
		      );



bool updateRTErrors(
		    bool onlyOverlap,
		    map<string, vector< vector< pair<remappedRead, remappedRead> > > >::iterator itf,
		    int NB_ROUNDS,
		    map<string, GFFtranscript> &mt,
		    map<string, string> &mg,
		    char * tab_cpt,
		    char * tabBases,
		    FILE *nf_candidates,
		    long *NB_RT_ERRORS,
		    long *NB_POS_COVERED_FOR_RT_ERRORS
		    );


/*
bool updateIndelErrors(
		    bool onlyOverlap,
		    map<string, vector< vector< pair<remappedRead, remappedRead> > > >::iterator itf,
		    int NB_ROUNDS,
		    map<string, GFFtranscript> &mt,
		    map<string, string> &mg,
		    char * tab_cpt,
		    char * tabBases,
		    FILE *nf_candidates,
		    long *NB_INSERTIONS,
		    long *NB_DELETIONS,
		    long *NB_POS_COVERED_FOR_INDELS
		    );

*/

void displayRTcandidate(map<string, vector< vector< pair<remappedRead, remappedRead> > > >::iterator it, int nbRounds, FILE *nf, char c1, char c2);

bool allMembersMapOnSameTranscript(vector< pair<remappedRead, remappedRead> > &v);



// FUNCTION TO DEAL WITH SEQUENCE QUALITY SCORES:



int updateTranscriptionErrors(
			     vector< vector< pair<remappedRead, remappedRead> > > &vFam,
			     double maxProbaError,
			     int nbMaxSoftEdits,
			     mapOfCoverage &mcov,
			     long minCov,
			     double minProp,
			     map<string, string> &mg,
			     map<string, GFFtranscript> &mt,
			     char *tabBases,
			     char *tab_cpt,
			     FILE *nf_candidates, 
			     long *NB_TR_ERRORS, 
			     long *NB_POS_COVERED_FOR_TR_ERRORS,
			     map<string, vector<long> > &mObs,
			     map<char, pair<long,long> > &mObsByBase
			     );


void getHardAndSoftConsensus(
			     vector< vector< pair<remappedRead, remappedRead> > > &vFam, 
			     map<long, pair<string, string> > &mExtendedConsensus,
			     map<long, pair<char, double> > &mhc, 
			     map<long, char> &msc,
			     pair<remappedRead, remappedRead> &prr,
			     map<long, pair< pair<char,char> , pair<char,char> > > &mi,
			     char *tabBases,
			     char *tab_cpt
			     );


void getExtendedConsensus(
			  vector< vector< pair<remappedRead, remappedRead> > > &vFam, 
			  map<long, pair<string, string> > &mc, 
			  pair<remappedRead, remappedRead> &prr,
			  map<long, pair< pair<char,char> , pair<char,char> > > &mi,
			  char *tab_cpt);



double getProbaHCisFalse(string &excQual);

void getIntersection(remappedRead &read, remappedRead &mate, map<long, pair< pair<char,char> , pair<char,char> > > &mi, char *tab_cpt);


void displayCandidate(
		      FILE *nf_candidates, 
		      vector< vector< pair<remappedRead, remappedRead> > > &vFam, 
		      char gc, 
		      char cc, 
		      string extendedConsensus,
		      string excQual,
		      double probaError,
		      long cov, 
		      double prop, 
		      pair<remappedRead, remappedRead> &prr,
		      long gPos, 
		      long tPos);




int updateRTerrors(
		   vector< vector< pair<remappedRead,remappedRead> > > &vFam,
		   double maxProbaError,
		   char *tabBases,
		   char *tab_cpt,
		   long *NB_RT_ERRORS,
		   long *NB_POS_COVERED_FOR_RT,
		   long *NB_NOT_REF_BASE,
		   long tabObsForRT[4],
		   long tabRterrors[4][4],
		   FILE *nf_candidates,
		   map<string,string> &mg
		   );


///////////////////////////////////////////////////////////////////////////////////////
// Retired version of the function
// int updateRTerrors(
//		 vector< vector< pair<remappedRead, remappedRead> > > &vFam,
//		 double maxProbaError,
//		 map<string, string> &mg,
//		 char *tabBases,
//		 char *tab_cpt,
//		 long *NB_RT_ERRORS,
//		 long *NB_POS_COVERED_FOR_RT,
//		 FILE *nf_candidates
//		 );
///////////////////////////////////////////////////////////////////////////////////////

void displayRTcandidate(
			FILE *nf_candidates,
			vector< vector< pair<remappedRead, remappedRead> > > &vFam,
			char c1,
			char c2,
			long gPos,
			pair<remappedRead,remappedRead> &pp,
			double p1,
			double p2
			);


void displayMap(FILE *nf, map<long, pair<string, string> > &m);
void displayListOfPosition(
			   map<string, vector< pair<long, long> > > &mapOfAllPos, 
			   map<string, vector< pair<remappedRead, remappedRead> > > * mr,
			   int nbRounds);

void displayAllIDs(
		   map<string, vector< pair<remappedRead, remappedRead> > > *mNbBarcodes, 
		   string barcode,
		   bool *roundHasPos, 
		   int nbRounds, 
		   FILE *nf);


void init_mObs(map<string, vector<long> > &mObs, map<string, string> &mg);
void display_mObs(map<string, vector<long> > &mObs, FILE *nf);

void sum_mapObs(map<string, vector<long> >&mt, map<string, vector<long> > *m, int MAX_SIZE);
void load_mapObs(map<string, vector<long> >&m, char *fname, bool header);

void init_mCodons(map<string, vector<long> > &m, char *tab_bases);

void display_mCodons(FILE *nf, map<string, vector<long> > &m);

#endif
