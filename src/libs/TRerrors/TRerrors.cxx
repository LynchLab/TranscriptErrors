#ifndef _TRerrors_CXX_
#define _TRerrors_CXX_

#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include <algorithm>
#include <vector>

#include "TRerrors.hxx"



long updateTranscriptionErrors(
			       long minCov,
			       double minProp,
			       bool onlyOverlap, 
			       map<string, vector< vector< pair<remappedRead, remappedRead> > > >::iterator itf, 
			       int nbRounds, 
			       map<string, GFFtranscript> &mt,
			       map<string, string> &mg, 
			       char *tab_cpt, 
			       mapOfCoverage &mcov,
			       char *tabBases,
			       FILE *nf_candidates, 
			       long *NB_TR_ERRORS, 
			       long *NB_POS_COVERED_FOR_TR_ERRORS,
			       map<string, vector<long> > &mObs
			       ){


  long res = 0;

  map<long, char> mhc;
  map<long, char> msc;

  vector< pair<remappedRead, remappedRead> > vFamily;

  int sizeFamily;
  int iRound;
  int withinRoundFamilySize;
  vector<int> vRoundSize;

  for(iRound=0, sizeFamily=0 ; iRound<nbRounds ; iRound++){

    withinRoundFamilySize = (itf->second)[iRound].size();

    // DEBOG :
    //if(iRound>0){
    //fprintf(stderr, "\t");
    //}
    //fprintf(stderr, "{");

    if( withinRoundFamilySize > 0 ){
      vRoundSize.push_back(withinRoundFamilySize);
      for(int ii=0 ; ii<withinRoundFamilySize ; ii++){

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// DEBOG !!!!!
	// !!!!!! TO BE REMOVEd LATTER !!!!!!!!!
	/*
	if( (itf->second)[iRound].at(ii).first.vDeletions_.size() > 0 || (itf->second)[iRound].at(ii).second.vDeletions_.size() > 0 ){
	  return 0;
	}
	if(ii>0){
	  fprintf(stderr, " , ");
	}
	fprintf(stderr, "%s", (itf->second)[iRound].at(ii).first.qName_.c_str());
	*/
	// !!!! END DEBOG !!!!!!!!!!!!!!!!!!!!!

	vFamily.push_back((itf->second)[iRound].at(ii) );
      }
      sizeFamily++;
    }
    //fprintf(stderr, "}");
  }

  //fprintf(stderr, "%s\n", vFamily[0].first.qName_.c_str());
  //fprintf(stderr, "\n");
  //fflush(stderr);

  int minRoundSize = vRoundSize[0];
  for(int ii=0 ; ii<vRoundSize.size() ; ii++){
    if(vRoundSize[ii] < minRoundSize){
      minRoundSize = vRoundSize[ii];
    }
  }

  onlyOverlap = true;
  if(vFamily.size()>3 || minRoundSize>1 || sizeFamily==3){
    onlyOverlap = false;
  }


  //onlyOverlap = true; // DEBOG

  getConsensus(vFamily, mhc, msc, tabBases, tab_cpt, onlyOverlap);

  //debogConsensus(vFamily, mhc, msc);
  
  res = updateTranscriptionErrors(
				  onlyOverlap,
				  minCov, 
				  minProp,
				  mhc, 
				  msc,
				  itf, 
				  mt, 
				  mg, 
				  tab_cpt, 
				  mcov,
				  tabBases,
				  nf_candidates, 
				  NB_TR_ERRORS, 
				  NB_POS_COVERED_FOR_TR_ERRORS,
				  mObs
				  );

  return res;

}


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
			       map<string, vector<long> > & mObs
			       ){


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Just a little bit of code to know what transcript / chromosome I'm looking at
  int nbRounds = (itf->second).size();
  int ii;
  for(ii=0 ; ii < nbRounds && (itf->second)[ii].size()==0 ; ii++)
    ;

  if(ii==nbRounds){
    fprintf(stderr, "ERROR: REACHED END OF VECTOR !!!\n");
    exit(1);
  }

  remappedRead rtmp = (itf->second)[ii].at(0).first;
  remappedRead mtmp = (itf->second)[ii].at(0).second;
  string chrID = rtmp.rName_;

  bool cpt = false;
  if( rtmp.strand_ != "+"){
    cpt = true;
  }

  string tID = rtmp.transcriptID_;
  // END OF THE FIRST PART USED TO FIND ON WHICH TRANSCRIPT / CHROMOSOME I'M CURRENTLY WORKING
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  long gPos;
  char gc;


  // I make sure that the alignment from the soft consensus is perfect (all bases are the same as the reference genome)
  //bool softAlnIsPerfect = true;
  int nbMismatchesInSoftAln = 0;
  map<long, char>::iterator itsc;
  for(itsc=msc.begin() ; itsc!=msc.end() ; itsc++){
    gPos = itsc->first;
    gc = mg[chrID][gPos];
    if(cpt==true){ gc = tab_cpt[gc]; }
    if( gc != itsc->second ){
      //softAlnIsPerfect = false;
      //return;
      nbMismatchesInSoftAln++;
    }
  }


  if(nbMismatchesInSoftAln>1){
    return 0;
  }


  int posInMap = 0;
  int mapSize = mhc.size();

  map<long, char>::iterator itc;
  for(itc=mhc.begin() ; itc!=mhc.end() ; itc++){
    char c = itc->second;
    posInMap++;
    if( c != 'N' ){

      gPos = itc->first;

      bool inRead = true;
      long readCoord = rtmp.genomeToReadCoords(gPos);
      if(readCoord == -1){
	inRead = false;
	readCoord = mtmp.genomeToReadCoords(gPos);
      }

      if( onlyOverlap==true ||  (posInMap>2 && posInMap<(mapSize-2) && readCoord>4 && (readCoord<86 || (inRead==false && readCoord<96) ) ) ){
	gc = mg[chrID][gPos];
	if(cpt==true){ gc = tab_cpt[gc]; }

	long cov = 0;
	double prop = 0.0;
	
	if(minCov>0 || minProp>0.0){
	  mcov.getCovAndProp(chrID, gPos, cpt, gc, tabBases, &cov, &prop);
	}

	if( cov >= minCov && prop >= minProp ){
	
	  (*NB_POS_COVERED_FOR_TR_ERRORS)++;
	  mObs[chrID].at(gPos)++;
	
	  if( gc != c ){
	    GFFtranscript transcript = mt[tID];
	    long tPos = transcript.genomicToSplicedTranscript(gPos);
	    displayCandidate(nf_candidates, gc, c,cov, prop, chrID, gPos, tID, tPos, itf, posInMap, mapSize);
	    fprintf(nf_candidates, "read coordinates :  %d\n", readCoord);
	    fflush(nf_candidates);
	    (*NB_TR_ERRORS)++;
	  }
	}// End of if coverage is ok
      }// End of if pos in map > 1
    }
  }
  return 1;
}


void getExtendedConsensus(vector< pair<remappedRead, remappedRead> > &vr, map<long, string> &mc, char *tabBases, char *tab_cpt, bool onlyOverlap){

  vector< pair<remappedRead, remappedRead> >::iterator itr;

  for(itr=vr.begin() ; itr!=vr.end() ; itr++){
    map<long, pair<char, char> > mi;
    getIntersection( itr->first, itr->second, mi, tab_cpt); 
    map<long, pair<char, char> >::iterator iti;
    for(iti=mi.begin() ; iti!=mi.end() ; iti++){
      long pos = iti->first;
      pair<char, char> pp = iti->second;

      if(onlyOverlap==true){ // If I look only at overlapping region, I include both read and mate base (so that I include a N in the consensus string sequence if I'm not in the overlapping region
	mc[pos].push_back(pp.first);
	mc[pos].push_back(pp.second);
      } else {
	if(pp.first!='N'){
	  mc[pos].push_back(pp.first);
	}
	if(pp.second!='N'){
	  mc[pos].push_back(pp.second);
	}
      }
    }
  }
}



void getCollapsedConsensus(vector< pair<remappedRead, remappedRead> > &vr, map<long, string> &mc, char *tabBases, char *tab_cpt){

  vector< pair<remappedRead, remappedRead> >::iterator itr;

  for(itr=vr.begin() ; itr!=vr.end() ; itr++){
    map<long, pair<char, char> > mi;
    getIntersection( itr->first, itr->second, mi, tab_cpt); 
    map<long, pair<char, char> >::iterator iti;
    for(iti=mi.begin() ; iti!=mi.end() ; iti++){
      long pos = iti->first;
      pair<char, char> pp = iti->second;

      if(pp.first==pp.second){
	mc[pos].push_back(pp.first);
      } else {
	mc[pos].push_back('N');
      }
    }
  }
}


void getConsensus(vector< vector< pair<remappedRead, remappedRead> > > &vFam, map<long, char> &mhc, map<long, char> &msc, char *tabBases, char *tab_cpt, bool onlyOverlap){

  vector< pair<remappedRead, remappedRead> > vtmp;

  int NB_ROUND = vFam.size();
  for(int iRound=0 ; iRound<NB_ROUND ; iRound++){
    int roundSize = vFam[iRound].size();
    for(int j=0 ; j<roundSize ; j++){
      vtmp.push_back(vFam[iRound].at(j));
    }
  }

  getConsensus(vtmp, mhc, msc, tabBases, tab_cpt, onlyOverlap);

}


void getConsensus(vector< pair<remappedRead, remappedRead> > &vr, map<long, char> &mhc, map<long, char> &msc, char *tabBases, char *tab_cpt, bool onlyOverlap){

  map<long, string> mc;
  getExtendedConsensus(vr, mc, tabBases, tab_cpt, onlyOverlap);

  map<long, string>::iterator it;
  for(it=mc.begin() ; it!=mc.end() ; it++){
    long pos = it->first;
    char hc, sc;
    getHardAndSoftConsensus(it->second, &hc, &sc, tabBases);
    mhc[pos] = hc;
    msc[pos] = sc;
  }
}




void getHardAndSoftConsensus(string &s, char *hc, char *sc, char *tabBases){

  int i, j;

  if( s.size()==0 ){
    (*hc) = 'N';
    (*sc) = 'N';
    return;
  }

  //if( s.find('-')!=string::npos ){
  //fprintf(stderr, "DEBOG : %s\n", s.c_str());
  //}

  int NB_BASES = 5;
  char tabBasesTmp[NB_BASES]; // = {'A', 'T', 'C', 'G', '-'};
  for(i=0 ; i<4 ; i++){
    tabBasesTmp[i] = tabBases[i];
  }
  tabBasesTmp[i] = '-';

  int lg = s.size();
  int tc[NB_BASES]; // = {0, 0, 0, 0};
  for(i=0 ; i<NB_BASES ; i++){
    tc[i] = 0;
  }

  for(i=0 ; i<lg ; i++){
    char base = s[i];
    for(j=0 ; j<NB_BASES && tabBasesTmp[j]!=base ; j++)
      ;
    if(j<NB_BASES){
      tc[j]++;
    }
  }

  int imax, max;
  for(i=0, max=-1, imax=-1 ; i<NB_BASES ; i++){
    if(tc[i]>max){
      max = tc[i];
      imax = i;
    }
  }

  if(max==lg){
    (*hc) = tabBasesTmp[imax];
  } else {
    (*hc) = 'N';
  }

  if(imax>=0){
    (*sc) = tabBasesTmp[imax];
  } else {
    (*sc) = 'N';
  }

}


int getIntersection(remappedRead &read, remappedRead &mate, map<long, pair<char, char> > &mi, char *tab_cpt){

  int i, posOnRead;
  long start, end;

  int nbOverlap = 0;
  long offset = 1;

  string readSeqWithIndels = read.getSeqWithIndels();
  string mateSeqWithIndels = mate.getSeqWithIndels();

  if( readSeqWithIndels=="" || mateSeqWithIndels=="" ){
    //fprintf(stderr, "ERROR in computing seqWithIndel for %s\n", read.qName_.c_str());
    //fprintf(stderr, "readSeqWithIndels: %s\n", readSeqWithIndels.c_str());
    //fprintf(stderr, "mateSeqWithIndels: %s\n", mateSeqWithIndels.c_str());
    return 0;
  }


  bool cpt = false;
  if(read.strand_=="-"){ 
    cpt = true;
    offset = -1;
  }


  //map<long, pair<char, char> > m;

  int nbSegRead = read.vmaps_.size();
  int nbSegMate = mate.vmaps_.size();

  posOnRead = 0;
  for(i=0 ; i<nbSegRead ; i++){
    long start = read.vmaps_[i].first;
    long lg = read.vmaps_[i].second;
    long end = (start + lg);
    if(cpt==true){ end = start - lg; }

    for(long pos=start ; pos!=end ; pos+=offset){

      char gc = readSeqWithIndels[posOnRead];

      pair<char, char> pp;
      pp.first = gc;
      pp.second = 'N';

      mi[pos] = pp;

      posOnRead++;
    }
  }


  posOnRead = 0;
  for(i=0 ; i<nbSegMate ; i++){
    long start = mate.vmaps_[i].first;
    long lg = mate.vmaps_[i].second;
    long end = (start + lg);
    if(cpt==true){ end = start - lg; }

    for(long pos=start ; pos!=end ; pos+=offset){
      char gc = mateSeqWithIndels[posOnRead];

      if( mi.find(pos)!=mi.end() ){
	mi[pos].second = gc;
	nbOverlap++;
      } else {
	pair<char, char> pp;
	pp.first = 'N';
	pp.second = gc;
	mi[pos] = pp;
      }

      posOnRead++;
    }
  }

  return nbOverlap;  
}


void prepareConsensus(map<long, char> &mc1, map<long, char> &mc){

  map<long, char>::iterator it;
  for(it=mc1.begin() ; it!=mc1.end() ; it++){
    long pos = it->first;
    char c = it->second;
    mc[pos] = c;
  }

}


void addConsensus(map<long, char> &mcs, map<long, char> &mc){

  map<long, char>::iterator it;
  for(it=mc.begin() ; it!=mc.end() ; it++){
    long pos = it->first;
    char c = it->second;
    map<long, char>::iterator it_tmp = mcs.find(pos);
    if( it_tmp == mcs.end() || it_tmp->second != c ){
      mc[pos] = 'N';
    }
  }

}


/*
void updateMapFamilies(
		       int iRound,
			map<string, vector< vector< pair< remappedRead, remappedRead> > > > & mf,
			map<string, vector< pair<remappedRead, remappedRead> > > & mBarcodesAtPosition,
			vector<string> & vBarcodesSeen,
			int nbRounds
			){

  map<string, vector< pair<remappedRead, remappedRead> > >::iterator it;
  for(it=mBarcodesAtPosition.begin() ; it!=mBarcodesAtPosition.end() ; it++){
    string barcode = it->first;
    if( find(vBarcodesSeen.begin(), vBarcodesSeen.end(), barcode)==vBarcodesSeen.end() ){
      for(int i=0 ; i<nbRounds ; i++){
	vector< pair<remappedRead, remappedRead> > vtmp;
	mf[barcode].push_back(vtmp);
      }
      vBarcodesSeen.push_back(barcode);
    }
    
    int nbRedondantPairs = (it->second).size();
    for(int ii=0 ; ii<nbRedondantPairs ; ii++){
      mf[barcode][iRound].push_back( (it->second)[ii] );
    }    
    
  }
  
}
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
		      char c){

  int iRound;

  fprintf(nf, 
	  "%s\t%d\t%c\t%c",
	  chrID.c_str(),
	  pos,
	  gc,
	  c
	  );

  for(iRound=0 ; iRound<nbRounds ; iRound++){
    if( roundHasPos[iRound]==true ){
      fprintf(nf, "\t%s", mNbBarcodes[iRound][barcode][0].first.qName_.c_str());
    }
  }

  fprintf(nf, "\n");

}
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
		      ){

  fprintf(
	  nf_candidates,
	  "%c\t%c\t%d\t%.4f\t%s\t%d\t%s\t%d\tposInMap=%d\tmapSize=%d",
	  gc,
	  c,
	  cov,
	  prop,
	  chrID.c_str(),
	  gPos,
	  tID.c_str(),
	  tPos,
	  posInMap,
	  mapSize
	  );

  int nbRounds = (itf->second).size();
  for(int ii=0 ; ii<nbRounds ; ii++){
    fprintf(nf_candidates, "\t{");
    int nbTmp = (itf->second)[ii].size();
    for(int jj=0 ; jj<nbTmp ; jj++){
      fprintf(nf_candidates, "%s", (itf->second)[ii][jj].first.qName_.c_str());
      if(jj<(nbTmp-1)){
	fprintf(nf_candidates, " , ");
      }
    }
    fprintf(nf_candidates, "}");
    //if( (itf->second)[ii].size() > 0 ){
    //fprintf(nf_candidates, "%s", (itf->second)[ii][0].first.qName_.c_str());
    //}
  }

  fprintf(nf_candidates, "\n");

}



/*

void getConsensus(vector< map<long, char> > &vmc, map<long, char> &mc){

  int nbSeq = vmc.size();

  prepareConsensus(vmc[0], mc);

  for(int i=1 ; i<nbSeq ; i++){
    addConsensus(vmc[i], mc);
  }
}


void getConsensus(vector< pair<remappedRead, remappedRead> > &vpr, map<long, char> &mc, char *tab_cpt, bool onlyOverlap){

  int nbr = vpr.size();
  vector< map<long, char> > vmc;
  for(int i=0 ; i<nbr ; i++){
    map<long, char> mc_tmp;
    getConsensus(vpr[i].first, vpr[i].second, mc_tmp, tab_cpt, onlyOverlap);
    vmc.push_back(mc_tmp);
  }

  getConsensus(vmc, mc);

}


void getConsensus(remappedRead &read, remappedRead &mate, map<long, char> &mc, char *tab_cpt, bool onlyOverlap){

  map<long, pair<char, char> > mi;
  getIntersection(read, mate, mi, tab_cpt);
  getConsensus(mi, mc, onlyOverlap);

}


void getConsensus(map<long, pair<char, char> > &mi, map<long, char> &mc, bool onlyOverlap){
  map<long, pair<char, char> >::iterator it;
  for(it=mi.begin() ; it!=mi.end() ; it++){
    long pos = it->first;
    pair<char, char> pp = it->second;

    if( onlyOverlap == true ){
      if(pp.first==pp.second){
	mc[pos] = pp.first;
      } else {
	mc[pos] = 'N';
      }
    } else {
      if(pp.first==pp.second && pp.first!='N'){
	mc[pos] = pp.first;
      } else {
	if(pp.first=='N'){
	  mc[pos] = pp.second;
	}
	if(pp.second=='N'){
	  mc[pos] = pp.first;
	}
	if(pp.first!='N' && pp.second!='N'){
	  mc[pos] = 'N';
	}
      }
    }
  }
}


*/

/*
//deprecated
bool updateRTErrors(
		    bool onlyOverlap,
		    map<string, vector< vector< pair<remappedRead, remappedRead> > > >::iterator itf,
		    int nbRounds,
		    map<string, GFFtranscript> &mt,
		    map<string, string> &mg,
		    char * tab_cpt,
		    char * tabBases,
		    FILE *nf_candidates,
		    long *NB_RT_ERRORS,
		    long *NB_POS_COVERED_FOR_RT_ERRORS
		    ){


  vector< pair<remappedRead, remappedRead> > vp_tmp;

  for(int iRound = 0 ; iRound<nbRounds  ; iRound++){
    int withinRoundFamilySize = (itf->second)[iRound].size();
    if( withinRoundFamilySize > 0 ){
      for(int ii=0 ; ii<withinRoundFamilySize ; ii++){
	vp_tmp.push_back( (itf->second)[iRound].at(ii) );
      }
    }
  }

  if(vp_tmp.size() != 2){
    return false;
  }

  if( allMembersMapOnSameTranscript(vp_tmp) == false ){
    return false;
  }

  map<long, string> mcs;
  getCollapsedConsensus(vp_tmp, mcs, tabBases, tab_cpt);

  int posInOverlap;

  map<long, string>::iterator it;
  for(it=mcs.begin(), posInOverlap=0 ; it!=mcs.end() ; it++){
    string s = it->second;
    if( s[0]!='N' && s[1]!='N' ){

      if( posInOverlap > 2 ){
	(*NB_POS_COVERED_FOR_RT_ERRORS)++;
	if( s[0]==s[1] ){
	  ;
	} else {
	  (*NB_RT_ERRORS)++;
	  displayRTcandidate(itf, nbRounds, nf_candidates, s[0], s[1]);
	}
      }
      posInOverlap++;
    }
  }

  return true;
}
*/

void displayRTcandidate(map<string, vector< vector< pair<remappedRead, remappedRead> > > >::iterator itf, int nbRounds, FILE *nf, char c1, char c2){

  fprintf(nf, "%c\t%c", c1, c2);
  
  for(int iRound = 0 ; iRound<nbRounds  ; iRound++){
    int withinRoundFamilySize = (itf->second)[iRound].size();
    if( withinRoundFamilySize > 0 ){
      fprintf(nf, "\t%s", (itf->second)[iRound][0].first.qName_.c_str());
    } else {
      fprintf(nf, "\tNO_READ");
    }
  }


  fprintf(nf, "\n");

}


bool allMembersMapOnSameTranscript(vector< pair<remappedRead, remappedRead> > &v){

  string tID = v[0].first.transcriptID_;

  int nb = v.size();

  for(int i=0 ; i<nb ; i++){
    if( v[i].first.transcriptID_ != tID || v[i].second.transcriptID_ != tID ){
      return false;
    }
  }
  return true;
}



/*
bool updateIndelErrors(
		    bool onlyOverlap,
		    map<string, vector< vector< pair<remappedRead, remappedRead> > > >::iterator itf,
		    int nbRounds,
		    map<string, GFFtranscript> &mt,
		    map<string, string> &mg,
		    char * tab_cpt,
		    char * tabBases,
		    FILE *nf_candidates,
		    long *NB_INSERTIONS,
		    long *NB_DELETIONS,
		    long *NB_POS_COVERED_FOR_INDELS
		       ){


  vector< pair<remappedRead, remappedRead> > vtmp_family;

  bool familyContainsIndels = false;

  for(int iRound = 0 ; iRound<nbRounds  ; iRound++){
    int withinRoundFamilySize = (itf->second)[iRound].size();

    if(withinRoundFamilySize>1){ return false; }

    if( withinRoundFamilySize > 0 ){
      for(int ii=0 ; ii<withinRoundFamilySize ; ii++){
	vtmp_family.push_back( (itf->second)[iRound].at(ii) );

	// !!! IN THIS VERSION, I ALLOW ONLY FAMILIES WITH ONE INDEL
	if( (itf->second)[iRound].at(ii).first.vInsertions_.size() == 1 ||
	    (itf->second)[iRound].at(ii).first.vDeletions_.size() == 1 ||
	    (itf->second)[iRound].at(ii).second.vInsertions_.size() == 1 ||
	    (itf->second)[iRound].at(ii).second.vDeletions_.size() == 1 ){

	  familyContainsIndels = true;

	}
      }
    }
  }

  map<long, pair<char, char> > mi;
  getIntersection(vtmp_family[0].first, vtmp_family[1].second, mi, tab_cpt);

  if(familyContainsIndels == false){
    (*NB_POS_COVERED_FOR_INDELS) += mi.size();
    return false;
  }

  int familySize = vtmp_family.size();

  //map< pair<long, long> , vector< pair<int, int> > > mDel;
  //map< pair<long, long> , vector< pair<int, int> > > mIns;

  vector< pair<long, long> > vIns;
  vector< pair<long, long> > vDel;

  int nbi = vtmp_family[0].first.vInsertions_.size();
  int posOfIns = 0; // 1=first ; 2=second ; 3=overlapping region

  if(nbi>1){ return false; }

  if( nbi==1 ){ 
    posOfIns = 1;
  }

  if( nbi==0 ){
    nbi = vtmp_family[0].second.vInsertions_.size();
    posOfIns = 2;
  }

  if( nbi==1 ){

    pair<long, long> pi;
    if(posOfIns==1){ 
      pi = vtmp_family[0].first.vInsertions_[0];
      int posInRead = vtmp_family[0].first.genomeToReadCoords(pi.first);
      if(posInRead<15 || posInRead>85){
	return false;
      }
    } else {
      pi = vtmp_family[0].second.vInsertions_[0];
      int posInMate = vtmp_family[0].second.genomeToReadCoords(pi.first);
      if(posInMate<15 || posInMate>85){
	return false;
      }
    }


    if( vtmp_family[0].first.coversPosition(pi.first) && vtmp_family[0].second.coversPosition(pi.first) ){

      posOfIns = 3;

      if( vtmp_family[0].first.asInsertion(pi)==false || vtmp_family[0].second.asInsertion(pi)==false ){
	return false;
      }

    }

    for(int ff=0 ; ff<familySize ; ff++){

      if( posOfIns==1 ){
	if( vtmp_family[ff].first.asInsertion(pi)==false || vtmp_family[ff].first.vInsertions_.size()>1 || vtmp_family[ff].second.vInsertions_.size()>0 ){
	  return false;
	}
      }

      if( posOfIns==2 ){
	if( vtmp_family[ff].second.asInsertion(pi)==false || vtmp_family[ff].second.vInsertions_.size()>1 || vtmp_family[ff].first.vInsertions_.size()>0 ){
	  return false;
	}
      }

      if( posOfIns==3 ){
	if( vtmp_family[ff].first.asInsertion(pi)==false || vtmp_family[ff].first.vInsertions_.size()>1 || 
	    vtmp_family[ff].second.asInsertion(pi)==false || vtmp_family[ff].second.vInsertions_.size()>1 ){
	  return false;
	}
      }
    }

    vIns.push_back(pi);
  }



  int nbd = vtmp_family[0].first.vDeletions_.size();
  int posOfDel = 0; // 1=first ; 2=second ; 3=overlapping region

  if(nbd>1){ return false; }

  if( nbd==1 ){ 
    posOfDel = 1;
  }

  if( nbd==0 ){
    nbd = vtmp_family[0].second.vDeletions_.size();
    posOfDel = 2;
  }

  if( nbd==1 ){

    pair<long, long> pd;
    if(posOfDel==1){ 
      pd = vtmp_family[0].first.vDeletions_[0];
      int posInRead = vtmp_family[0].first.genomeToReadCoords(pd.first);
      if(posInRead<15 || posInRead>85){
	return false;
      }
    } else {
      pd = vtmp_family[0].second.vDeletions_[0];
      int posInMate = vtmp_family[0].second.genomeToReadCoords(pd.first);
      if(posInMate<15 || posInMate>85){
	return false;
      }
    }


    if( vtmp_family[0].first.coversPosition(pd.first) && vtmp_family[0].second.coversPosition(pd.first) ){

      posOfDel = 3;

      if( vtmp_family[0].first.asDeletion(pd)==false || vtmp_family[0].second.asDeletion(pd)==false ){
	return false;
      }

    }

    for(int ff=0 ; ff<familySize ; ff++){

      if( posOfDel==1 ){
	if( vtmp_family[ff].first.asDeletion(pd)==false || vtmp_family[ff].first.vDeletions_.size()>1 || vtmp_family[ff].second.vDeletions_.size()>0 ){
	  return false;
	}
      }

      if( posOfDel==2 ){
	if( vtmp_family[ff].second.asDeletion(pd)==false || vtmp_family[ff].second.vDeletions_.size()>1 || vtmp_family[ff].first.vDeletions_.size()>0 ){
	  return false;
	}
      }

      if( posOfDel==3 ){
	if( vtmp_family[ff].first.asDeletion(pd)==false || vtmp_family[ff].first.vDeletions_.size()>1 || 
	    vtmp_family[ff].second.asDeletion(pd)==false || vtmp_family[ff].second.vDeletions_.size()>1 ){
	  return false;
	}
      }
    }

    vDel.push_back(pd);
  }



  //(*NB_POS_COVERED_FOR_INDELS) += mi.size();

  if( vIns.size() == 1 ){

    (*NB_INSERTIONS)++;

    remappedRead read_tmp = vtmp_family[0].first;
    pair<long, long> pi = vIns[0];

    fprintf(nf_candidates, 
	    "I\t%s\t%d\t%d\t%s",
	    read_tmp.rName_.c_str(),
	    pi.first,
	    pi.second,
	    read_tmp.transcriptID_.c_str()
	    );

    for(int iRound = 0 ; iRound<nbRounds  ; iRound++){
      vector< pair<remappedRead, remappedRead> > vtmp = (itf->second)[iRound];
      if( vtmp.size() == 0 ){
	fprintf(nf_candidates, "\tNO_READ");
      } else {
	fprintf(nf_candidates, "\t%s", vtmp[0].first.qName_.c_str());
      }
    }

    fprintf(nf_candidates, "\n");
    
  }


  if( vDel.size() == 1 ){

    (*NB_DELETIONS)++;

    remappedRead read_tmp = vtmp_family[0].first;
    pair<long, long> pd = vDel[0];

    fprintf(nf_candidates, 
	    "D\t%s\t%d\t%d\t%s",
	    read_tmp.rName_.c_str(),
	    pd.first,
	    pd.second,
	    read_tmp.transcriptID_.c_str()
	    );

    for(int iRound = 0 ; iRound<nbRounds  ; iRound++){
      vector< pair<remappedRead, remappedRead> > vtmp = (itf->second)[iRound];
      if( vtmp.size() == 0 ){
	fprintf(nf_candidates, "\tNO_READ");
      } else {
	fprintf(nf_candidates, "\t%s", vtmp[0].first.qName_.c_str());
      }
    }

    fprintf(nf_candidates, "\n");
    
  }


  return true;
}
*/


		       
void debogConsensus(vector< pair<remappedRead, remappedRead> > &vFamily, map<long, char> &mhc, map<long, char> &msc){

  fprintf(stderr, 
	  "%s\n%s\n%s\n",
	  vFamily[0].first.qName_.c_str(),
	  vFamily[0].first.seq_.c_str(),
	  vFamily[0].second.seq_.c_str()
	  );

  fprintf(stderr, "Hard consensus : \n");
  map<long, char>::iterator it;
  for(it=mhc.begin() ; it!=mhc.end() ; it++){
    fprintf(stderr, 
	    "%d\t%c\n", 
	    it->first,
	    it->second
	    );
  }

  fprintf(stderr, "\nSoft consensus : \n");
  for(it=msc.begin() ; it!=msc.end() ; it++){
    fprintf(stderr, 
	    "%d\t%c\n", 
	    it->first,
	    it->second
	    );
  }


  fprintf(stderr, "-------------------------------\n");

}





// VERSION WITH TAKING INTO ACCOUNT THE SEQUENCING QUALITY SCORES:
// ---------------------------------------------------------------


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
			     map<char, pair<long, long> > &mObsByBase
			     ){


  int familySize = 0;

  int MIN_POS_READ = 5;
  int MAX_POS_READ = 103;
  int MIN_POS_MATE = 0;
  int MAX_POS_MATE = 95;

  int totalSize = 0;
  int nbRound = vFam.size();
  for(int iRound=0 ; iRound<nbRound ; iRound++){
    totalSize += vFam[iRound].size();
    if(vFam[iRound].size()>0){
      familySize++;
    }
  }

  // !! DEBOG !!
  //if(totalSize>6){
  //return 0;
  //}

  map<long, pair<char, double> > mhc;
  map<long, char> msc;

  pair<remappedRead, remappedRead> prr;
  map<long, pair< pair<char,char> , pair<char,char> > > mi;

  map<long, pair<string, string> > mExtendedConsensus;

  getHardAndSoftConsensus(vFam, mExtendedConsensus, mhc, msc, prr, mi, tabBases, tab_cpt);

  int readMapSize = prr.first.vmaps_.size();
  int mateMapSize = prr.second.vmaps_.size();
  int overSeq;

  if( prr.first.strand_ == "+" ){
    overSeq = prr.first.vmaps_[0].first - prr.second.vmaps_[0].first;
  } else {
    overSeq =  prr.second.vmaps_[0].first - prr.first.vmaps_[0].first;

    /*
    if(readMapSize>1 && mateMapSize>1){
      fprintf(stderr, "Read map:\n");
      for(int ii=0 ; ii<readMapSize ; ii++){
	fprintf(stderr, "%d : %d\n", prr.first.vmaps_[ii].first, prr.first.vmaps_[ii].second);
      }
      fprintf(stderr, "Mate map:\n");
      for(int ii=0 ; ii<mateMapSize ; ii++){
	fprintf(stderr, "%d : %d\n", prr.second.vmaps_[ii].first, prr.second.vmaps_[ii].second);
      }
      exit(1);
    }
    */
  }

  if( overSeq > 0 ){
    MIN_POS_MATE = overSeq + 4;
    fprintf(stderr, "Mate goes past read (overSeq=%d) : %s\n", overSeq, prr.first.qName_.c_str());
    //return 0;
  }



  string chrID = prr.first.rName_;
  string tID = prr.first.transcriptID_;
  bool cpt = false;
  if( prr.first.strand_ == "-" ){ cpt = true; };


  int nbSoftEdits = 0;
  map<long, char>::iterator itsc;
  for(itsc=msc.begin() ; itsc!=msc.end() ; itsc++){

    long gPos = itsc->first;

    long readCoord = prr.first.genomeToReadCoords(gPos);
    long mateCoord = prr.second.genomeToReadCoords(gPos);

    if( (readCoord<0 && mateCoord<0) || ( readCoord>=0 && (readCoord<MIN_POS_READ || readCoord>MAX_POS_READ) ) || (mateCoord>=0 && (mateCoord<MIN_POS_MATE || mateCoord>MAX_POS_MATE) ) ){
      continue;
    }

    char gc = mg[chrID].at(gPos);
    if( cpt==true ){ gc = tab_cpt[gc]; }

    if(gc!=itsc->second){ nbSoftEdits++; }

  }

  if( nbSoftEdits>nbMaxSoftEdits ){
    return 0;
  }

  bool mito = false;
  if(prr.first.rName_=="MtDNA"){
    mito = true;
  }

  int posInMap;
  int mapHCsize = mhc.size();
  map<long, pair<char, double> >::iterator ithc;
  for(ithc = mhc.begin(), posInMap=0 ; ithc != mhc.end() ; ithc++, posInMap++){

    if( posInMap>2 && posInMap<(mapHCsize-2) ){

      long gPos = ithc->first;
      string extendedConsensus = mExtendedConsensus[gPos].first;
      string excQual = mExtendedConsensus[gPos].second;


      //bool inRead = true;
      long readCoord = prr.first.genomeToReadCoords(gPos);
      long mateCoord = prr.second.genomeToReadCoords(gPos);

      //if( readCoord>0 && mateCoord<0){
      //fprintf(stderr, "before if: %d - %d\n", readCoord, mateCoord);
      //fflush(stderr);
      //}

      if( (readCoord<0 && mateCoord<0) || ( readCoord>=0 && (readCoord<MIN_POS_READ || readCoord>MAX_POS_READ) ) || (mateCoord>=0 && (mateCoord<MIN_POS_MATE || mateCoord>MAX_POS_MATE) ) ){
	//continue;
	;
      } else {

	double probaError = (ithc->second).second;

	if( probaError <= maxProbaError ){

	  char cc = (ithc->second).first;
	  char gc = mg[chrID].at(gPos);
	  if( cpt==true ){
	    gc = tab_cpt[gc];
	  }

	  long cov = 0;
	  double prop = 0.0;
	
	  if(minCov>0 || minProp>0.0){
	    mcov.getCovAndProp(chrID, gPos, cpt, gc, tabBases, &cov, &prop);
	  }

	  if( cov >= minCov && prop >= minProp ){
	
	    NB_POS_COVERED_FOR_TR_ERRORS[familySize]++;
	    mObs[chrID].at(gPos)++;

	    if(mito==false){
	      mObsByBase[gc].first++;
	    } else {
	      mObsByBase[gc].second++;
	    }

	    //if( readCoord>0 && mateCoord<0){
	    //fprintf(stderr, "%d - %d\n", readCoord, mateCoord);
	    //fflush(stderr);
	    //}

	    if( gc != cc ){
	      GFFtranscript transcript = mt[tID];
	      long tPos = transcript.genomicToSplicedTranscript(gPos);
	      displayCandidate(nf_candidates, vFam, gc, cc, extendedConsensus, excQual, probaError, cov, prop, prr, gPos, tPos);
	      //fprintf(nf_candidates, "prr: %s %d (readCoord:%d , mateCoord:%d)\n", prr.first.strand_.c_str(), overSeq, readCoord , mateCoord);
	      fflush(nf_candidates);
	      NB_TR_ERRORS[familySize]++;
	    }// End of if gc!=cc (= end of if I see a transcription error)
	    
	  }// End of if the coverage/prop is enough to look at this position
	
	}// End of if the position in read/mate is ok

      }// End of if the probability of erroneous base call is low enough

    }// End of if posInMap is ok

  }// End of the loop that goes through all the positions in the hard consensus map ( for(ithc = mhc.begin() ; ithc != mhc.end() ; ithc++){ ... )

  return 1;

}


void getHardAndSoftConsensus(
			     vector< vector< pair<remappedRead, remappedRead> > > &vFam, 
			     map<long, pair<string, string> > &mExtendedConsensus,
			     map<long, pair<char, double> > &mhc, 
			     map<long, char> &msc, 
			     pair<remappedRead, remappedRead> &prr,
			     map<long, pair< pair<char,char> , pair<char,char> > > &mi,
			     char *tabBases,
			     char *tab_cpt
			     ){

  //map<long, pair<string, string> > mExtendedConsensus;
  getExtendedConsensus(vFam, mExtendedConsensus, prr, mi, tab_cpt);

  map<long, pair<string, string> >::iterator it;
  for(it=mExtendedConsensus.begin() ; it!=mExtendedConsensus.end() ; it++){
    long gPos = it->first;
    string excSeq = (it->second).first;
    string excQual = (it->second).second;

    char sc, hc;
    getHardAndSoftConsensus(excSeq, &hc, &sc, tabBases);

    msc[gPos] = sc;

    pair<char, double> pp;
    pp.first = hc;
    if( hc=='N' ){
      pp.second = 1.0;
    } else {
      pp.second = getProbaHCisFalse(excQual);
    }

    mhc[gPos] = pp;

  }

}


double getProbaHCisFalse(string &excQual){

  int iTotalQual = 0;
  int nb = excQual.size();

  double proba = 3.0;

  for(int i=0 ; i<nb ; i++){
    char c = excQual[i];
    int q = (int)c - 33;
    iTotalQual += q;
    proba = proba * ( ( (double)1.0/(double)3.0 ) * (double)1/pow(10, q/(double)10.0 ) );
  }

  //double proba = (double)1/pow(10,((double)iTotalQual/(double)10.0) );


  return proba;
}


void getExtendedConsensus(
			  vector< vector< pair<remappedRead, remappedRead> > > &vFam, 
			  map<long, pair<string, string> > &mc, 
			  pair<remappedRead, remappedRead> &prr,
			  map<long, pair< pair<char,char> , pair<char,char> > > &mi,
			  char *tab_cpt){


  bool initialized = false;

  int nbRounds = vFam.size();

  for(int iRound=0 ; iRound<nbRounds ; iRound++){
    int roundSize = vFam[iRound].size();
    for(int jj=0 ; jj<roundSize ; jj++){
      //pair<remappedRead, remappedRead> prrTmp = vFam[iRound].at(jj);
      prr = vFam[iRound].at(jj);
      //map<long, pair< pair<char, char>, pair<char, char> > > mi;
      getIntersection(prr.first, prr.second, mi, tab_cpt);

      map<long, pair< pair<char,char> , pair<char,char> > >::iterator it;
      for(it=mi.begin() ; it!=mi.end() ; it++){
	long gPos = it->first;

	char c1 = (it->second).first.first;
	char c2 = (it->second).second.first;

	char q1 = (it->second).first.second;
	char q2 = (it->second).second.second;

	if( c1=='N' && c2=='N' )
	  continue;

	if( c1=='N' ){
	  mc[gPos].first.push_back(c2);
	  mc[gPos].second.push_back(q2);
	  continue;
	}

	if( c2=='N' ){
	  mc[gPos].first.push_back(c1);
	  mc[gPos].second.push_back(q1);
	  continue;
	}

	mc[gPos].first.push_back(c1);
	mc[gPos].first.push_back(c2);
	mc[gPos].second.push_back(q1);
	mc[gPos].second.push_back(q2);

      }

    }
  }

}





void getIntersection(remappedRead &read, remappedRead &mate, map<long, pair< pair<char,char> , pair<char,char> > > &mi, char *tab_cpt){

  int i, posOnRead;
  long start, end;

  long offset = 1;

  string readSeqWithIndels = read.getSeqWithIndels(false);
  string mateSeqWithIndels = mate.getSeqWithIndels(false);

  string readQualWithIndels = read.getSeqWithIndels(true);
  string mateQualWithIndels = mate.getSeqWithIndels(true);


  if( readSeqWithIndels=="" || mateSeqWithIndels=="" ){
    //fprintf(stderr, "ERROR in computing seqWithIndel for %s\n", read.qName_.c_str());
    //fprintf(stderr, "readSeqWithIndels: %s\n", readSeqWithIndels.c_str());
    //fprintf(stderr, "mateSeqWithIndels: %s\n", mateSeqWithIndels.c_str());
    return;
  }


  bool cpt = false;
  if(read.strand_=="-"){ 
    cpt = true;
    offset = -1;
  }

  int nbSegRead = read.vmaps_.size();
  int nbSegMate = mate.vmaps_.size();

  posOnRead = 0;
  for(i=0 ; i<nbSegRead ; i++){
    long start = read.vmaps_[i].first;
    long lg = read.vmaps_[i].second;
    long end = (start + lg);
    if(cpt==true){ end = start - lg; }

    for(long pos=start ; pos!=end ; pos+=offset){

      char gc = readSeqWithIndels[posOnRead];
      char q = readQualWithIndels[posOnRead];

      pair< pair<char,char> , pair<char,char> > pp;

      pp.first.first = gc;
      pp.first.second = q;

      pp.second.first = 'N';
      pp.second.second = '#';

      mi[pos] = pp;

      posOnRead++;
    }
  }


  posOnRead = 0;
  for(i=0 ; i<nbSegMate ; i++){
    long start = mate.vmaps_[i].first;
    long lg = mate.vmaps_[i].second;
    long end = (start + lg);
    if(cpt==true){ end = start - lg; }

    for(long pos=start ; pos!=end ; pos+=offset){

      char gc = mateSeqWithIndels[posOnRead];
      char q = mateQualWithIndels[posOnRead];

      if( mi.find(pos)!=mi.end() ){
	mi[pos].second.first = gc;
	mi[pos].second.second = q;
      } else {
	pair< pair<char,char> , pair<char,char> > pp;
	pp.first.first = 'N';
	pp.first.second = '#';

	pp.second.first = gc;
	pp.second.second = q;

	mi[pos] = pp;
      }

      posOnRead++;
    }
  }
  
}








void init(map<long, pair<string, int> >&mc, map<long, pair< pair<char,int> , pair<char,int> > > &mi){

  map<long, pair< pair<char,int> , pair<char,int> > >::iterator it;
  for(it=mi.begin() ; it!=mi.end() ; it++){
    long gPos = it->first;
    pair<string, int> pp;
    pp.first = "";
    pp.second = 0;
    mc[gPos] = pp;
  }

}



void displayCandidate(
		      FILE *nf_candidates, 
		      vector< vector< pair<remappedRead, remappedRead> > > &vFam, 
		      char gc, 
		      char cc, 
		      string extendedConsensus,
		      string excQual,
		      double proba,
		      long cov, 
		      double prop, 
		      pair<remappedRead, remappedRead> &prr,
		      long gPos, 
		      long tPos){

  string chrID = prr.first.rName_;
  string tID = prr.first.transcriptID_;
  string strand = prr.first.strand_;

  fprintf(nf_candidates, "%c\t%c", gc, cc);

  //fprintf(nf_candidates, "\t%s\t%s\t%.6f", extendedConsensus.c_str(), excQual.c_str(), (double)1000*proba);

  fprintf(nf_candidates, "\t%s\t%d\t%s\t%s", chrID.c_str(), gPos, tID.c_str(), strand.c_str());

  fprintf(nf_candidates, "\tcov:%d\tprop:%.6f", cov, prop);

  int nbRounds = vFam.size();
  for(int iRound = 0 ; iRound<nbRounds ; iRound++){
    int roundSize = vFam[iRound].size();
    fprintf(nf_candidates, "\t{");
    for(int jj=0 ; jj<roundSize ; jj++){
      if(jj>0){ 
	fprintf(nf_candidates, " , ");
      }
      fprintf(nf_candidates, "%s", vFam[iRound].at(jj).first.qName_.c_str());
    }
    fprintf(nf_candidates, "}");
  }

  long readPos = prr.first.genomeToReadCoords(gPos);
  long matePos = prr.second.genomeToReadCoords(gPos);
  fprintf(nf_candidates, "\t%d\t%d", readPos, matePos);

  fprintf(nf_candidates, "\n");

  //fprintf(stderr, "%c->%c (proba: %.8f)\n", gc, cc, proba);

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////
// FROM HERE: FUNCTIONS TO ESTIMATE THE RT ERROR RATE:

int updateRTerrors(
		   vector< vector< pair<remappedRead,remappedRead> > > &vFam,
		   double maxProbaError,
		   char *tabBases,
		   char *tab_cpt,
		   long *NB_RT_ERRORS,
		   long *NB_POS_COVERED_FOR_RT,
		   long *NB_NOT_REF_BASE,
		   long tabObsForRT[4],
		   long tabRTerrors[4][4],
		   FILE *nf_candidates,
		   map<string,string> &mg
		   ){

  int i, iRound, nbRounds;
  int i_refBase, i_alternateBase;
  
  map<long, pair< pair<char,char> , pair<char,char> > > mi;
  map<long, pair<string, string> > mExtendedConsensus;
  pair<remappedRead, remappedRead> prr;
  map<long, pair<char, double> > mhc;
  map<long, char> msc;

  vector< map<long, pair<char, double> > > v_mhc;

  vector< vector< pair<remappedRead, remappedRead> > > vtmp;


  // The first thing to do is to make a vector of hard consensus for each round of RT. The result is stored in a vector(v_mhc)
  int nbInformativeRounds = 0;
  nbRounds = vFam.size();
  for(iRound = 0 ; iRound<nbRounds ; iRound++){
    if( vFam[iRound].size() > 0 ){
      vtmp.push_back(vFam[iRound]);
      getHardAndSoftConsensus(vtmp, mExtendedConsensus, mhc, msc, prr, mi, tabBases, tab_cpt);
      v_mhc.push_back(mhc);

      //fprintf(stderr, "Extended concensus:\n");
      //displayMap(stderr, mExtendedConsensus);

      mExtendedConsensus.clear();
      mhc.clear();
      msc.clear();
      mi.clear();

      nbInformativeRounds++;

      vtmp.clear();

    }
  }


  map<long, pair<char,double> >::iterator ithc;
  int mapSize = v_mhc[0].size();
  int posInMap = 0;

  for(ithc=v_mhc[0].begin(), posInMap=0 ; ithc!=v_mhc[0].end() ; ithc++, posInMap++){

    if(posInMap<1 || posInMap>(mapSize-2))
      continue;

    long gPos = ithc->first;

    //fprintf(stderr, "%ld -> %c\t%c\n", gPos, ithc->second.first, mg[prr.first.rName_].at(gPos) );
    //fprintf(stderr, "%ld -> %c", gPos, mg[prr.first.rName_].at(gPos) );

    int nbObsAtPos[5] = {0,0,0,0};
    double maxProba = ithc->second.second;

    for(iRound=0 ; iRound<nbInformativeRounds ; iRound++){
      pair<char,double> pp_tmp = v_mhc[iRound][gPos];
      char c = pp_tmp.first;
      double proba = pp_tmp.second;

      if(proba>maxProba){
	maxProba = proba;
      }

      for(i=0 ; i<4 && c!=tabBases[i] ; i++)
	;
      nbObsAtPos[i]++; // If c=='N' (or any other non A/T/C/G letter, then i==4)

    }

    if( maxProba>maxProbaError || nbObsAtPos[4]!=0 )
      continue;

    char refBase = mg[prr.first.rName_].at(gPos);
    if( prr.first.strand_=="-" ){
      refBase = tab_cpt[refBase];
    }

    int nbUniqueBases = 0;
    bool refBaseIsPresent = false;

    for(i=0 ; i<4 ; i++){
      if( nbObsAtPos[i]!=0 ){
	nbUniqueBases++;
	if( tabBases[i] == refBase ){
	  i_refBase = i;
	  refBaseIsPresent = true;
	} else {
	  i_alternateBase = i;
	}

      }
    }

    // If I have a family of size 3 and all 3 have a different base call, I skip the position (this situation should almost never happen)
    // nbUniqueBases<1 should never happen...
    if(nbUniqueBases<1 || nbUniqueBases>2)
      continue;

    if(refBaseIsPresent==false){
      (*NB_NOT_REF_BASE)++;
      continue;
    }

    (*NB_POS_COVERED_FOR_RT)++;
    tabObsForRT[i_refBase]++;

    if(nbUniqueBases==2){

      char alternateBase = tabBases[i_alternateBase];

      if(nbObsAtPos[i_alternateBase]!=1){
	fprintf(stderr, "PROBLEM!\n");
	displayRTcandidate(stderr, vFam, refBase, alternateBase, gPos, prr, 0.0, 0.0);
      }

      (*NB_RT_ERRORS)++;

      // These are the errors as they appear on the mRNA/read. To obtain the real RT errors, it is necessary to take the reverse complementary version of the bases.
      tabRTerrors[i_refBase][i_alternateBase]++;

      // TO DO: output the candidate in the nf_candidate file.
      displayRTcandidate(nf_candidates, vFam, refBase, alternateBase, gPos, prr, 0.0, 0.0);
      fflush(nf_candidates);

    }// End of if(nbUniqueBases==2)


  }// End of the for loop that goes position by position in the hard consensus (   for(ithc=v_mhc[0].begin() ; ithc!=v_mhc[0].end() ; ithc++) )

  return 1;

}


/*
// Retired version
int updateRTerrors(
		 vector< vector< pair<remappedRead, remappedRead> > > &vFam,
		 double maxProbaError,
		 map<string, string> &mg,
		 char *tabBases,
		 char *tab_cpt,
		 long *NB_RT_ERRORS,
		 long *NB_POS_COVERED_FOR_RT,
		 FILE *nf_candidates
		 ){

  vector< vector< pair<remappedRead, remappedRead> > > v1;
  vector< vector< pair<remappedRead, remappedRead> > > v2;

  int nbRounds = vFam.size();

  int iRound, nbDone;
  for(iRound=0, nbDone=0 ; iRound<nbRounds && nbDone<2 ; iRound++){
  //for(iRound=1, nbDone=0 ; iRound<nbRounds && nbDone<2 ; iRound++){
    if(iRound!=-1){
      if( vFam[iRound].size() > 0 ){
	if(nbDone==0){
	  v1.push_back(vFam[iRound]);
	} else {
	  v2.push_back(vFam[iRound]);
	}
	nbDone++;
      }
    }
  }


  map<long, pair< pair<char,char> , pair<char,char> > > mi1;
  map<long, pair< pair<char,char> , pair<char,char> > > mi2;
  map<long, pair<string, string> > mExtendedConsensus1;
  map<long, pair<string, string> > mExtendedConsensus2;
  pair<remappedRead, remappedRead> prr;
  map<long, pair<char, double> > mhc1;
  map<long, pair<char, double> > mhc2;
  map<long, char> msc1;
  map<long, char> msc2;


  getHardAndSoftConsensus(v1, mExtendedConsensus1, mhc1, msc1, prr, mi1, tabBases, tab_cpt);
  getHardAndSoftConsensus(v2, mExtendedConsensus2, mhc2, msc2, prr, mi2, tabBases, tab_cpt);


  map<long, pair<char,double> >::iterator ithc;
  for(ithc=mhc1.begin() ; ithc!=mhc1.end() ; ithc++){
    long gPos = ithc->first;
    char c1 = mhc1[gPos].first;
    double p1 = mhc1[gPos].second;
    char c2 = mhc2[gPos].first;
    double p2 = mhc2[gPos].second;

    //    if(c1!='N' && c2!='N'){
    if( (c1=='A' || c1=='T' || c1=='C' || c1=='G') && (c2=='A' || c2=='T' || c2=='C' || c2=='G') ){

      if(p1<maxProbaError && p2<maxProbaError){
	(*NB_POS_COVERED_FOR_RT)++;
	if(c1!=c2){
	  (*NB_RT_ERRORS)++;
	  displayRTcandidate(nf_candidates, vFam, c1, c2, p1, p2);
	  fflush(nf_candidates);
	}
      }
    }
  }

  return 1;
}
*/

void displayRTcandidate(
			FILE *nf_candidates,
			vector< vector< pair<remappedRead, remappedRead> > > &vFam,
			char c1,
			char c2,
			long gPos,
			pair<remappedRead, remappedRead> &pp,
			double p1,
			double p2
			){

  string chrID = pp.first.rName_;
  string tID = pp.first.transcriptID_;

  fprintf(nf_candidates, 
	  "%c\t%c\t%s\t%d\t%s\t%s\tcov:%.4f\tprop:%.4f",
	  c1,
	  c2,
	  chrID.c_str(),
	  gPos,
	  tID.c_str(),
	  pp.first.strand_.c_str(),
	  p1*(double)1000.0,
	  p2*(double)1000.0
	  );

  int nbRound = vFam.size();
  for(int iRound=0 ; iRound<nbRound ; iRound++){
    fprintf(nf_candidates, "\t{");
    int roundSize = vFam[iRound].size();
    for(int ii=0 ; ii<roundSize ; ii++){
      if(ii>0){
	fprintf(nf_candidates, " , ");
      }
      fprintf(nf_candidates, "%s", vFam[iRound].at(ii).first.qName_.c_str());
    }
    fprintf(nf_candidates, "}");
  }

  fprintf(nf_candidates, "\n");

}


void displayMap(FILE *nf, map<long, pair<string, string> > &m){

  map<long, pair<string, string> >::iterator it;
  for(it=m.begin() ; it!=m.end() ; it++){
    fprintf(
	    nf, 
	    "%d\t%s\t%s\n", 
	    it->first,
	    it->second.first.c_str(),
	    it->second.second.c_str()
	    );
  }

}


void display_mObs(map<string, vector<long> > &mObs, FILE *nf){

  vector< pair<string, long> > vObs;

  map<string, vector<long> >::iterator it;
  for(it=mObs.begin() ; it!=mObs.end() ; it++){
    string chrID = it->first;
    long lg = (it->second).size();

    long tot = 0;

    for(long ll=0 ; ll<lg ; ll++){
      long obs = (it->second).at(ll);
      fprintf(nf, 
	      "%s\t%d\t%d\n", 
	      chrID.c_str(),
	      ll,
	      obs
	      );

      tot += obs;

    }

    pair<string, long> pp;
    pp.first = chrID;
    pp.second = tot;
    vObs.push_back(pp);
  }

  int nb = vObs.size();
  for(int ii=0 ; ii<nb ; ii++){
    fprintf(nf, "#%s\t%d\n", vObs[ii].first.c_str(), vObs[ii].second);
  }

}


void displayListOfPosition(
			   map<string, vector< pair<long, long> > > &mapOfAllPos, 
			   map<string, vector< pair<remappedRead, remappedRead> > > * mr,
			   int nbRounds
			   ){


  int i, iRound;
  string transcriptID;

  map<string, vector< pair<long, long> > >::iterator it_all_pos;
  for(it_all_pos=mapOfAllPos.begin() ; it_all_pos!=mapOfAllPos.end() ; it_all_pos++){
    transcriptID = it_all_pos->first;
    int nbPos = (it_all_pos->second).size();

    map< pair<long, long> , vector< pair<remappedRead, remappedRead> > > mgr_tmp[nbRounds];
    bool roundHasTranscript[3];

    for(iRound=0 ; iRound<nbRounds ; iRound++){
      if( mr[iRound].find(transcriptID) != mr[iRound].end() ){
	getReadsGroupedByPosition(mr[iRound][transcriptID], mgr_tmp[iRound]);
	roundHasTranscript[iRound] = true;
      } else {
	roundHasTranscript[iRound] = false;
      }
    }

    for(i=0 ; i<nbPos ; i++){

      pair<long, long> pp_tmp = (it_all_pos->second)[i];
      int nbRoundsSharingPos = 0;

      for(iRound=0 ; iRound<nbRounds ; iRound++){
	if( roundHasTranscript[iRound]==true && mgr_tmp[iRound].find(pp_tmp)!=mgr_tmp[iRound].end() ){
	  nbRoundsSharingPos++;
	}
      }

      fprintf(stdout, 
	      "%s\t%d\t%d\t%d\n", 
	      transcriptID.c_str(), 
	      pp_tmp.first, 
	      pp_tmp.second,
	      nbRoundsSharingPos
	      );
    }
  }
}


void displayAllIDs(
		   map<string, vector< pair<remappedRead, remappedRead> > > *mNbBarcodes, 
		   string barcode,
		   bool *roundHasPos, 
		   int nbRounds, 
		   FILE *nf){


  for(int iRound=0 ; iRound<nbRounds ; iRound++){
    map<string, vector< pair<remappedRead, remappedRead> > > mb = mNbBarcodes[iRound];
    int nbp = mb[barcode].size();
    //for(int ip=0 ; ip<nbp ; ip++){
    //fprintf(stdout, "ROUND_%d\t%s\n", iRound, mb[barcode][ip].first.qName_.c_str());
    //}
    fprintf(stdout, "%s", mb[barcode][0].first.qName_.c_str());
    if(iRound<(nbRounds-1)){
      fprintf(stdout, "\t");
    }
  }

  fprintf(stdout, "\n");

}


void init_mObs(map<string, vector<long> > &mObs, map<string, string> &mg){

  map<string, string>::iterator it;

  for(it=mg.begin() ; it!=mg.end() ; it++){
    long lg = (it->second).size();
    string chrID = it->first;
    vector<long> vtmp(lg, 0);
    mObs[chrID] = vtmp;
  }

}


void sum_mapObs(map<string, vector<long> >&mt, map<string, vector<long> > *m, int MAX_SIZE){

  int i;

  map<string, vector<long> >::iterator it;
  for(it=m[0].begin() ; it!=m[0].end() ; it++){
    string chrID = it->first;
    long lg = it->second.size();
    for(long pos=0 ; pos<lg ; pos++){
      long total = 0;
      for(int ff=1 ; ff<=MAX_SIZE ; ff++){
	long val = m[ff][chrID].at(pos);
	total += val;
      }
      mt[chrID].at(pos) = total;
    }
  }

}

void load_mapObs(map<string, vector<long> >&m, char *fname, bool header){

  string line, chrID, stmp;
  long nbObs;
  long pos = 0;

  string previousChrID = "";

  ifstream fic(fname);
  if(header==true){
    getline(fic, line);
  }

  while( getline(fic, line) ){

    if( line[0] != '#' ){
      istringstream iss(line);
      getline(iss, chrID, '\t');

      getline(iss, stmp, '\t');
      pos = atol(stmp.c_str());

      getline(iss, stmp, '\t');
      nbObs = atol(stmp.c_str());

      m[chrID].at(pos) = nbObs;
    }

  }

}

void init_mCodons(map<string, vector<long> > &m, char *tab_bases){

  int i, j, k;
  char codon[4];
  codon[3] = '\0';
  
  for(i=0 ; i<4 ; i++){
    codon[0] = tab_bases[i];
    for(j=0 ; j<4 ; j++){
      codon[1] = tab_bases[j];
      for(k=0 ; k<4 ; k++){
	codon[2] = tab_bases[k];
	string sc(codon);
	vector<long> vv(3, 0);
	m.insert(map<string, vector<long> >::value_type(sc,vv));
      }
    }
  }
  

}

void display_mCodons(FILE *nf, map<string, vector<long> > &m){

  map<string, vector<long> >::iterator it;

  for(it=m.begin() ; it!=m.end() ; it++){
    string codon = it->first;
    vector<long> vv = it->second;
    int vSize = vv.size();
    fprintf(nf, "%s", codon.c_str());
    for(int ii=0 ; ii<vSize ; ii++){
      fprintf(nf, "\t%d", vv[ii]);
    }
    fprintf(nf, "\n");
  }

}

#endif

