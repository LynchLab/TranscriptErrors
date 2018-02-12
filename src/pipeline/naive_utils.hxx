#ifndef _NAIVE_UTILS_HXX
#define _NAIVE_UTILS_HXX

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>

#include "UTILS_FASTA/utils_fasta.hxx"
#include "CircSeqConsensus/CS_consensus.hxx"


bool fillVectorFromArgs(std::vector<std::string> &vArgs, std::string &sArgs, char sep);
int utils_findPosInMate(std::ifstream &ficBlatPos, std::string &id);
//int naive_utils_getDebreakedCs(CS_consensus &cs, std::string &id, std::ifstream &ficBlatPos, std::ifstream &ficCsPos, std::ifstream &ficReads, std::ifstream &ficMates, bool paired, char *tab_cpt);

int naive_utils_getDebreakedCs(
			       CS_consensus &cs, 
			       string &id, 
			       ifstream &ficBlatPos, string &lineBlat,
			       ifstream &ficCsPos,
			       ifstream &ficReads, 
			       ifstream &ficMates, 
			       bool paired, 
			       char *tab_cpt,
			       int *csLgParam,
			       int *alnLgParam,
			       int *breakPosParam
			       );



int naive_utils_getDebreakedCs(
			       CS_consensus &cs, 
			       string &id, 
			       ifstream &ficBlatPos,
			       ifstream &ficCsPos,
			       ifstream &ficReads, 
			       ifstream &ficMates, 
			       bool paired, 
			       char *tab_cpt
			       );

#endif
