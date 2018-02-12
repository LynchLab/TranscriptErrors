#ifndef _REMAPPEDREAD_UTILS_HXX_
#define _REMAPPEDREAD_UTILS_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>


#include "UTILS_FASTA/utils_fasta.hxx"
#include "remappedRead.hxx"
#include "TRerrors/TRerrors.hxx"


using namespace std;


void initMapGroupedPairedReadsBy_mRNA(
				      char *fname, 
				      map<string, vector< pair<remappedRead, remappedRead> > > & m, 
				      char minLetter, 
				      char maxLetter, 
				      string strand, 
				      bool includeSequence, 
				      bool computeSeqWithIndel, 
				      bool includeQuality);

void initMapGroupedPairedReadsBy_mRNA(
				      char *fname, 
				      map<string, vector< pair<remappedRead, remappedRead> > > & m, 
				      vector<string> & vIDs,
				      char minLetter, 
				      char maxLetter, 
				      string strand, 
				      bool includeSequence, 
				      bool computeSeqWithIndel, 
				      bool includeQuality);



void getReadsGroupedByPosition(vector< pair<remappedRead, remappedRead> > &v, map< pair<long, long> , vector< pair<remappedRead, remappedRead> > > &m);
void getReadsGroupedByPosition(vector< pair<remappedRead, remappedRead> > &v, map<long , vector< pair<remappedRead, remappedRead> > > &m);




void loadMapBarcodesAtPosition(
		       map< pair<long, long> , vector< pair<remappedRead, remappedRead> > >::iterator it, 
		       map< string, vector< pair<remappedRead, remappedRead> > > & mNbBarcodes,
		       map<string, string> & mBarcodes
		       );


void loadMapBarcodesAtPosition(
			       vector< pair<remappedRead, remappedRead> > &v,
			       map<string, vector< pair<remappedRead, remappedRead> > > & mNbBarcodes,
			       map<string, string> & mBarcodes
			       );


void loadBarcodes(char *fname, map<string, string> &mb);

long getMapOfAllPos(map<string, vector< pair<long, long> > > & mapOfAllPos, map<string, vector< pair<remappedRead, remappedRead> > > * mr, int nbRounds);
long getMapOfAllPos(map<string, vector<long> > & mapOfAllPos, map<string, vector< pair<remappedRead, remappedRead> > > * mr, int nbRounds);



void updateMapFamilies(
		       int iRound,
			map<string, vector< vector< pair< remappedRead, remappedRead> > > > & mf,
			map<string, vector< pair<remappedRead, remappedRead> > > & mBarcodesAtPosition,
			vector<string> & vBarcodesSeen,
			int nbRounds
		       );


void getFamiliesAtPosition(
			   int nbRounds, 
			   map<string, vector< vector< pair<remappedRead, remappedRead> > > > &mFamilies, 
			   map< pair<long, long>, vector< pair<remappedRead, remappedRead> > > * mgr_tmp, 
			   pair<long, long> &currentPos,  
			   map<string, string> *mBarcodes,
			   bool *roundHasTranscript
			   );


void getFamiliesAtPosition(
			   int nbRounds, 
			   map<string, vector< vector< pair<remappedRead, remappedRead> > > > &mFamilies, 
			   map< long, vector< pair<remappedRead, remappedRead> > > * mgr_tmp, 
			   long currentPos,  
			   map<string, string> *mBarcodes,
			   bool *roundHasTranscript
			   );




void displayFamily(FILE *nf, vector< vector< pair<remappedRead,remappedRead> > > &vFam);
void displayFamily( vector< vector< pair<remappedRead, remappedRead> > > &vFam, int NB_ROUND, FILE *nf);
void displayFamily(vector< vector< pair<remappedRead, remappedRead> > > &vFam, string barcode, int NB_ROUND, FILE *nfBarcodesInfo[], FILE *nfPairedInfo[]);


#endif
