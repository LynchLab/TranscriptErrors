#ifndef _REMAPPEDREAD_UTILS_CXX_
#define _REMAPPEDREAD_UTILS_CXX_
/*

remappedRead_utils.cxx

 */



#include "remappedRead_utils.hxx"



using namespace std;


void initMapGroupedPairedReadsBy_mRNA(
				      char *fname, 
				      std::map<std::string, std::vector< std::pair<remappedRead, remappedRead> > > & m, 
				      char minLetter, 
				      char maxLetter, 
				      string strand, 
				      bool includeSequence, 
				      bool computeSeqWithIndel, 
				      bool includeQuality
				      ){


  vector<string> vIDs;
  initMapGroupedPairedReadsBy_mRNA(fname, m, vIDs, minLetter, maxLetter, strand, includeSequence, computeSeqWithIndel, includeQuality);

}


void initMapGroupedPairedReadsBy_mRNA(
				      char *fname, 
				      std::map<std::string, std::vector< std::pair<remappedRead, remappedRead> > > & m, 
				      std::vector<std::string> & vIDs,
				      char minLetter, 
				      char maxLetter, 
				      string strand, 
				      bool includeSequence, 
				      bool computeSeqWithIndel, 
				      bool includeQuality
				      ){
  
  string line, s_read, s_mate;
  
  ifstream fic(fname);
  fprintf(stderr, "Grouping the reads by matching transcript ID ...");
  while( getline(fic, line) ){
    istringstream iss(line);
    getline(iss, s_read, '|');
    getline(iss, s_mate, '|');
    remappedRead read(s_read, includeSequence, computeSeqWithIndel, includeQuality);
    remappedRead mate(s_mate, includeSequence, computeSeqWithIndel, includeQuality);

    string idGeneRead = utils_fasta_extractGeneID(read.transcriptID_, '_');//getGeneID_fromTranscriptID(read.transcriptID_);;
    string idGeneMate = utils_fasta_extractGeneID(mate.transcriptID_, '_');//getGeneID_fromTranscriptID(mate.transcriptID_);;

    bool geneIsInList = false;

    if( vIDs.size()==0 ){
      if(idGeneRead[0]>=minLetter && idGeneRead[0]<=maxLetter && (strand=="b" || read.strand_==strand) && (strand=="b" || mate.strand_==strand) ){
	geneIsInList = true;
      }
    } else {
      if( find(vIDs.begin(), vIDs.end(), idGeneRead) != vIDs.end() ){
	geneIsInList = true;
      }
    }

    if(geneIsInList==true){
    
      if( idGeneRead != idGeneMate ){
	fprintf(stderr, "Warning, read and mate do not map to the same gene !\n%s-%s : %s-%s\n", read.qName_.c_str(), mate.qName_.c_str(), idGeneRead.c_str(), idGeneMate.c_str());
      }// else {  DEBOG !!!
      pair<remappedRead, remappedRead> p(read, mate);
      m[idGeneRead].push_back(p);
	//}
    }
  }
  fprintf(stderr, " done.\n");

}


void getReadsGroupedByPosition(vector< pair<remappedRead, remappedRead> > &v, map< pair<long, long> , vector< pair<remappedRead, remappedRead> > > &m){

  long nbr = v.size();
  for(long ll=0 ; ll<nbr ; ll++){
    long posRead = v[ll].first.vmaps_[0].first;
    long posMate = v[ll].second.vmaps_[0].first;
    pair<long, long> pp(posRead, posMate);
    m[pp].push_back(v[ll]);
  }
}



void getReadsGroupedByPosition(vector< pair<remappedRead, remappedRead> > &v, map<long , vector< pair<remappedRead, remappedRead> > > &m){

  long nbr = v.size();
  for(long ll=0 ; ll<nbr ; ll++){
    long posRead = v[ll].first.vmaps_[0].first;
    m[posRead].push_back(v[ll]);
  }

}


void loadMapBarcodesAtPosition(
		       map< pair<long, long> , vector< pair<remappedRead, remappedRead> > >::iterator it, 
		       map< string, vector< pair<remappedRead, remappedRead> > > & mNbBarcodes,
		       map<string, string> & mBarcodes
		       ){

  // it is an iterator in the map of positions
  int groupSize = (it->second).size();

  for(int vv=0 ; vv<groupSize ; vv++){
    pair<remappedRead, remappedRead> pp = (it->second)[vv];

    map<string,string>::iterator it = mBarcodes.find(pp.first.qName_);
    if(it==mBarcodes.end()){
      fprintf(stderr, "ERROR: NO BARCODE FOUND FOR '%s'\nABORTING.\n", pp.first.qName_.c_str());
      exit(4);
    }

    //string barcode = mBarcodes[pp.first.qName_];
    string barcode = it->second;
    mNbBarcodes[barcode].push_back(pp);
  }

}



void loadMapBarcodesAtPosition(
			       vector< pair<remappedRead, remappedRead> > &v,
			       map< string, vector< pair<remappedRead, remappedRead> > > & mNbBarcodes,
			       map<string, string> & mBarcodes
			       ){

  int groupSize = v.size();

  for(int ii=0 ; ii<groupSize ; ii++){
    pair<remappedRead, remappedRead> pp = v[ii];
    string barcode = mBarcodes[pp.first.qName_];
    mNbBarcodes[barcode].push_back(pp);
  }

}


void loadBarcodes(char *fname, map<string, string> &mb){

  string line, readID, barcode;

  ifstream ficBarcodes(fname);
  while( getline(ficBarcodes, line) ){
    istringstream iss(line);
    getline(iss, readID, '\t');
    getline(iss, barcode, '\t');
    mb.insert(map<string, string>::value_type(readID, barcode));
  }
  ficBarcodes.close();

}

long getMapOfAllPos(map<string, vector<long> > & mapOfAllPos, map<string, vector< pair<remappedRead, remappedRead> > > * mr, int nbRounds){

  long nbPosTotal = 0;
  int iRound;

  map<string, vector< pair<remappedRead, remappedRead> > >::iterator it_mr;
  for(iRound=0 ; iRound<nbRounds ; iRound++){
    for(it_mr=mr[iRound].begin() ; it_mr!=mr[iRound].end() ; it_mr++){
      string transcriptID = it_mr->first;
      map<string, vector< pair<remappedRead, remappedRead> > >::iterator it_mr_tmp = mr[iRound].find(transcriptID);
      if( it_mr_tmp != mr[iRound].end() ){
	map< long , vector< pair<remappedRead, remappedRead> > > mgr_tmp;
	getReadsGroupedByPosition(it_mr_tmp->second, mgr_tmp);
	
	map< long , vector< pair<remappedRead, remappedRead> > >::iterator it_tmp;
	for(it_tmp=mgr_tmp.begin() ; it_tmp!=mgr_tmp.end() ; it_tmp++){
	  long pos = it_tmp->first;
	  if( find(mapOfAllPos[transcriptID].begin(), mapOfAllPos[transcriptID].end(), pos) == mapOfAllPos[transcriptID].end() ){
	    mapOfAllPos[transcriptID].push_back(pos);
	    nbPosTotal++;
	  }
	}
      }
    }
  }// End ofthe loop used to record all the positions covered in the analysis.

  return nbPosTotal;
}


long getMapOfAllPos(map<string, vector< pair<long, long> > > & mapOfAllPos, map<string, vector< pair<remappedRead, remappedRead> > > * mr, int nbRounds){

  int iRound;
  long nbPosTotal = 0;

  map<string, vector< pair<remappedRead, remappedRead> > >::iterator it_mr;
  for(iRound=0 ; iRound<nbRounds ; iRound++){
    for(it_mr=mr[iRound].begin() ; it_mr!=mr[iRound].end() ; it_mr++){
      string transcriptID = it_mr->first;
      map<string, vector< pair<remappedRead, remappedRead> > >::iterator it_mr_tmp = mr[iRound].find(transcriptID);
      if( it_mr_tmp != mr[iRound].end() ){
	map< pair<long, long> , vector< pair<remappedRead, remappedRead> > > mgr_tmp;
	getReadsGroupedByPosition(it_mr_tmp->second, mgr_tmp);
	
	map< pair<long, long> , vector< pair<remappedRead, remappedRead> > >::iterator it_tmp;
	for(it_tmp=mgr_tmp.begin() ; it_tmp!=mgr_tmp.end() ; it_tmp++){
	  pair<long, long> pp_tmp = it_tmp->first;
	  if( find(mapOfAllPos[transcriptID].begin(), mapOfAllPos[transcriptID].end(), pp_tmp) == mapOfAllPos[transcriptID].end() ){
	    mapOfAllPos[transcriptID].push_back(pp_tmp);
	    nbPosTotal++;
	  }
	}
      }
    }
  }// End ofthe loop used to record all the positions covered in the analysis.

  return nbPosTotal;
}




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



void getFamiliesAtPosition(
			   int nbRounds, 
			   map<string, vector< vector< pair<remappedRead, remappedRead> > > > &mFamilies, 
			   map< pair<long, long> , vector< pair<remappedRead ,remappedRead> > > * mgr_tmp, 
			   pair<long, long> & currentPos, 
			   map<string, string> *mBarcodes, 
			   bool *roundHasTranscript){

  int iRound;

  vector<string> vBarcodesSeen;

  for(iRound=0 ; iRound<nbRounds ; iRound++){
    if( roundHasTranscript[iRound]==true && mgr_tmp[iRound].find(currentPos)!=mgr_tmp[iRound].end() ){
      map<string, vector< pair<remappedRead, remappedRead> > > mBarcodesAtPosition;
      loadMapBarcodesAtPosition(mgr_tmp[iRound].find(currentPos), mBarcodesAtPosition, mBarcodes[iRound]);      
      updateMapFamilies(iRound, mFamilies, mBarcodesAtPosition, vBarcodesSeen, nbRounds);
    }
  }
  
}

void getFamiliesAtPosition(
			   int nbRounds, 
			   map<string, vector< vector< pair<remappedRead, remappedRead> > > > &mFamilies, 
			   map< long , vector< pair<remappedRead ,remappedRead> > > * mgr_tmp, 
			   long currentPos, 
			   map<string, string> *mBarcodes, 
			   bool *roundHasTranscript){

  int iRound;

  vector<string> vBarcodesSeen;

  for(iRound=0 ; iRound<nbRounds ; iRound++){
    if( roundHasTranscript[iRound]==true && mgr_tmp[iRound].find(currentPos)!=mgr_tmp[iRound].end() ){
      map<string, vector< pair<remappedRead, remappedRead> > > mBarcodesAtPosition;
      loadMapBarcodesAtPosition(mgr_tmp[iRound][currentPos], mBarcodesAtPosition, mBarcodes[iRound]);
      updateMapFamilies(iRound, mFamilies, mBarcodesAtPosition, vBarcodesSeen, nbRounds);
    }
  }
  
}




void displayFamily(FILE *nf, vector< vector< pair<remappedRead,remappedRead> > > &vFam){

  int i, j;
  pair<remappedRead, remappedRead> prr;

  int nbRounds = vFam.size();
  int famSize = 0;

  for(i=0 ; i<nbRounds ; i++){
    if( vFam[i].size()>0 ){
      famSize++;
      prr = vFam[i].at(0);
    }
  }

  fprintf(nf, "%s", prr.first.transcriptID_.c_str());

  for(i=0 ; i<nbRounds ; i++){
    fprintf(nf, "\t{");
    int roundSize = vFam[i].size();
    for(j=0 ; j<roundSize ; j++){
      fprintf(nf, "%s", vFam[i].at(j).first.qName_.c_str());
    }
    fprintf(nf, "}");
  }
  fprintf(nf,"\n");

}

void displayFamily(vector< vector< pair<remappedRead, remappedRead> > > &vFam, string barcode, int NB_ROUND, FILE *nfBarcodesInfo[], FILE *nfPairedInfo[]){

  for(int iRound=0 ; iRound<NB_ROUND ; iRound++){
    int groupSize = vFam[iRound].size();
    for(int j=0 ; j<groupSize ; j++){
      // Printing out the paired reads
      vFam[iRound].at(j).first.fullPrint(nfPairedInfo[iRound], true);
      //fprintf(nfPairedInfo[iRound], "first");
      fprintf(nfPairedInfo[iRound], "|");
      vFam[iRound].at(j).second.fullPrint(nfPairedInfo[iRound], true);
      //fprintf(nfPairedInfo[iRound], "second");
      fprintf(nfPairedInfo[iRound], "\n");

      // Printing out the barcode association
      fprintf(nfBarcodesInfo[iRound], 
	      "%s\t%s\n", 
	      //"readID",
	      vFam[iRound].at(j).first.qName_.c_str(),
	      barcode.c_str()
	      );
    }
  }
}

void displayFamily( vector< vector< pair<remappedRead, remappedRead> > > &vFam, int NB_ROUND, FILE *nf){

  for(int i=0 ; i<NB_ROUND ; i++){
    if(i>0){
      fprintf(nf, "\t");
    }

    int groupSize = vFam[i].size();
    fprintf(nf, "{");
    for(int j=0 ; j<groupSize ; j++){
      if(j>0){
	fprintf(nf, " , ");
      }
      fprintf(nf, "%s", vFam[i].at(j).first.qName_.c_str());
    }

    fprintf(nf, "}");
  }
  fprintf(nf, "\n");
}






#endif

