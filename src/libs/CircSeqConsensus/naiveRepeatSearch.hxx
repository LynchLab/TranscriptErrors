#ifndef _NAIVE_REPEAT_SEARCH_HXX_
#define _NAIVE_REPEAT_SEARCH_HXX_

#include <iostream>
#include <fstream>
#include <sstream>


bool findConsensus(const std::string &consensus, const std::string &read, const std::string qual, const int maxMM);

int findFirstRepeatNaive(const std::string &s, const int minRepeatSize, const int maxMM);
int findSeqWithMM(const std::string &seqToFind, const std::string &seq, const int maxMM);

int findFirstOccurrenceNaive(const std::string &sToFind, const std::string &sToSearch, const int maxMM, int *nbMM, const int startSearchFromPos=0);


#endif
