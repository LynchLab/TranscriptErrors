#ifndef _GFFTRANSCRIPT_UTILS_HXX_
#define _GFFTRANSCRIPT_UTILS_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <map>
#include <vector>


#include "GFFtranscripts.hxx"

using namespace std;


void getAllTranscripts( 
		       std::map<std::string, GFFtranscript> & mst,  
		       std::map<std::string, std::vector<GFFentry> > & mvGFF, 
		       string sTAG,
		       int shift,
		       std::map<std::string, std::string> & mg, 
		       char *tab_cpt
			);


void getAllTranscripts(char *fname, int nbHeaderLines, map<string, GFFtranscript> &mst, string format, string sTAG, int shift, map<string, string> & mg, char *tab_cpt, vector<string> &vSources, vector<string> &vFeatures, vector<string> &vSeqNames);
void getAllTranscripts(char *fname, int nbHeaderLines, map<string, GFFtranscript> &mst, string format, string sTAG, int shift, map<string, string> & mg, char *tab_cpt);
void getAllTranscripts(char *fname, int nbHeaderLines, map<string, GFFtranscript> &mst, string format, string sTAG, int shift, vector<string> &vSources, vector<string> &vFeatures, vector<string> &vSeqNames);
void getAllTranscripts(char *fname, int nbHeaderLines, map<string, GFFtranscript> &mst, string format, string sTAG, int shift);

void getPosOfIntronsGenomicScale(map<string, vector<bool> > &mpi,  map<string, GFFtranscript> &mst, map<string, string> &mg, bool exonic);

void transcript2genomicCoords(GFFtranscript & gt, vector<long> & vg);

void addStartStop(std::map<std::string, GFFtranscript> &mst, std::vector<GFFentry> &vStarts, std::vector<GFFentry> &vStops, std::string sTag="transcript_id");

void orderTranscriptsByChromosomeAndStartPosition(map<string, GFFtranscript> &mt, map<string, vector<GFFtranscript> > &mvt);

void getAllTranscriptsAtPosition(vector<GFFtranscript> &vtr, map<string, GFFtranscript> &mt, string & chrID, long pos);
//void getAllTranscriptsAtPosition(vector<GFFtranscript> &vtr, map<string, vector<GFFtranscript> > &mvt, string chrID, long pos);// Not implemented yet
void getAllTranscriptsAtPosition(vector<GFFtranscript> &vtr, vector<GFFtranscript> &vt, long pos);


#endif
