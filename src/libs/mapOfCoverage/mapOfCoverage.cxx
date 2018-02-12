#ifndef _mapOfCoverage_CXX_
#define _mapOfCoverage_CXX_

#include "mapOfCoverage.hxx"



mapOfCoverage::mapOfCoverage(void){
  ;
}


void mapOfCoverage::init(map<string, string> &mg, bool verbose){

  map<string, string>::iterator it;
  for(it=mg.begin() ; it!=mg.end() ; it++){
    string chrID = it->first;
    long lg = (it->second).size();
    //vector<coverage> vtmp(lg, coverage((long)0));
    vector<coverage> vtmp;
    coverage cv( (long)0 );
    for(long ll=0 ; ll<lg ; ll++){
      vtmp.push_back(cv);
    }

    this->coverage_[chrID] = vtmp;
    if( verbose==true ){
      fprintf(stderr, " %s,", chrID.c_str());
    }
  }
}


void mapOfCoverage::updateFromPreviousMap(char *fname){


  string line;
  string chrID = "";
  string previousChrID = "";

  //string sval[NB_BASES_MOC];
  string sval_tmp;
  long tval[NB_BASES_MOC];

  long pos = 0;

  ifstream fic(fname);
  getline(fic, line); //Read the header line
  while( getline(fic, line) ){
    istringstream iss(line);
    getline(iss, chrID, '\t');

    if( chrID != previousChrID ){
      pos = 0;
      previousChrID = chrID;
    }

    // This part could is commented because should never happen (speeds up the execution time).
    //if( this->coverage_.find(chrID) == this->coverage_.end() ){
    //fprintf(stderr, "ERROR: '%s' not found in map of coverage. Interrupting the programm...\n", chrID.c_str());
    //exit(1);
    //}

    for(int ii=0 ; ii<NB_BASES_MOC ; ii++){
      getline(iss, sval_tmp, '\t');
      tval[ii] = atol(sval_tmp.c_str());
    }
    this->coverage_[chrID].at(pos).add(tval);

    pos++;
  }

}


void mapOfCoverage::updateFromRemapped(char *fname, map<string, string> &mg, char *tab_cpt, map<string, string> &mBarcodes, char *tabBases, bool verbose){
  map<string, vector< pair<remappedRead, remappedRead> > > m;
  bool includeSequence = true;
  bool computeSeqWithIndels = false;
  bool includeQuality = false;
  initMapGroupedPairedReadsBy_mRNA(fname, m, '0', 'z', "b", includeSequence, computeSeqWithIndels, includeQuality);
  this->updateFromRemapped(m, mg, tab_cpt, mBarcodes, tabBases, verbose);
  m.clear();
}


void mapOfCoverage::updateFromRemapped(map<string, vector< pair<remappedRead, remappedRead> > > & mr, map<string, string> &mg, char *tab_cpt, map<string, string> &mBarcodes, char *tabBases, bool verbose){

  long nbFamiliesDone = 0;

  map<string, vector< pair<remappedRead, remappedRead> > >::iterator it_mr;

  for(it_mr=mr.begin() ; it_mr!=mr.end() ; it_mr++){
    map< pair<long, long> , vector< pair<remappedRead, remappedRead> > > mrp;
    getReadsGroupedByPosition(it_mr->second, mrp);

    map< pair<long, long>, vector< pair<remappedRead, remappedRead> > >::iterator itp;
    for(itp=mrp.begin() ; itp!=mrp.end() ; itp++){
      map<string, vector< pair<remappedRead, remappedRead> > > mBarcodesAtPosition;
      loadMapBarcodesAtPosition(itp, mBarcodesAtPosition, mBarcodes);

      map<string, vector< pair<remappedRead, remappedRead> > >::iterator itb;
      for(itb=mBarcodesAtPosition.begin() ; itb!=mBarcodesAtPosition.end() ; itb++){
	map<long, char> mhc;
	map<long, char> msc;
	getConsensus(itb->second, mhc, msc, tabBases, tab_cpt, false);
	string chrID = (itb->second)[0].first.rName_;
	bool cpt = false;
	if( (itb->second)[0].first.strand_ == "-" ){
	  cpt = true;
	}

	this->update(msc, chrID, tab_cpt, cpt, tabBases);
	nbFamiliesDone++;

	if(verbose==true && nbFamiliesDone%100000 == 0 ){
	  fprintf(stderr, "%d families done so far ...\n", nbFamiliesDone);
	}

      }

    }

  }

}


void mapOfCoverage::update(map<long, char> &mc, string &chrID, char *tab_cpt, bool cpt, char *tabBases){

  map<string, vector<coverage> >::iterator itmcov_tmp = this->coverage_.find(chrID);

  if( itmcov_tmp == this->coverage_.end() ){
    fprintf(stderr, "ERROR: read.rName_ ('%s') not found. Aborting.\n", chrID.c_str());
    exit(2);
  }

  map<long, char>::iterator itmcons_tmp;
  for(itmcons_tmp=mc.begin() ; itmcons_tmp!=mc.end() ; itmcons_tmp++){
    long pos = itmcons_tmp->first;
    char c = itmcons_tmp->second;
    //if( cpt==true ){
    //c = tab_cpt[c];
    //}
    (itmcov_tmp->second)[pos].add(c, tabBases);
  }
}




void mapOfCoverage::updateFromFamily(vector< pair<remappedRead, remappedRead> > &vf, char *tab_cpt, char *tabBases){

  map<long, char> mhc;
  map<long, char> msc;
  getConsensus(vf, mhc, msc, tabBases, tab_cpt, false);

  string chrID = vf[0].first.rName_;

  bool cpt = false;
  if( vf[0].first.strand_ == "-" ){
    cpt = true;
  }

  this->update(msc, chrID, tab_cpt, cpt, tabBases);

}









void getCoverageAt(string chrID, long pos, coverage &c){
  ;
}

/////////////////////////////////////////////////////////////////////////////////////////

void mapOfCoverage::displayHeader(FILE *nf, char *tabBases){

  fprintf(nf, "chrID");
  for(int ii=0 ; ii<NB_BASES_MOC ; ii++){
    fprintf(nf, "\t%c", tabBases[ii]);
  }
  fprintf(nf, "\n");
}


void mapOfCoverage::display(FILE *nf, char *tabBases){

  this->displayHeader(nf, tabBases);

  map<string, vector<coverage> >::iterator it;
  for(it=this->coverage_.begin() ; it!=this->coverage_.end() ; it++){
    string chrID = it->first;
    long lg = (it->second).size();
    for(long ll=0 ; ll<lg ; ll++){
      coverage cc = (it->second)[ll];
      fprintf(nf, "%s", chrID.c_str());
      cc.display(nf);
      fprintf(nf, "\n");
    }
  }

}


void mapOfCoverage::display(FILE *nf, char *tabBases, map<string, string> &mg){

  //this->displayHeader(nf, tabBases);
  fprintf(nf, "%s\t%s\t%s", "chrID", "pos", "ref");

  for(int ii=0 ; ii<NB_BASES_MOC ; ii++){
    fprintf(nf, "\t%c", tabBases[ii]);
  }
  fprintf(nf, "\n");

  map<string, vector<coverage> >::iterator it;
  for(it=this->coverage_.begin() ; it!=this->coverage_.end() ; it++){
    string chrID = it->first;
    long lg = (it->second).size();

    for(long ll=0 ; ll<lg ; ll++){
      coverage cc = (it->second)[ll];
      fprintf(nf, "%s\t%d\t%c", chrID.c_str(), ll, mg[chrID].at(ll));
      cc.display(nf);
      fprintf(nf, "\n");
    }
  }

}


void mapOfCoverage::getCovAndProp(string &chrID, long gPos, bool cpt, char gc, char *tabBases, long *totalCov, double *prop){

  coverage cov = this->coverage_[chrID].at(gPos);

  long tmpTotal = 0;
  long covEq = 0;
  for(int ii=0 ; ii<NB_BASES_MOC ; ii++){
    long ltmp = cov.adjustedCov_[ii];
    if( tabBases[ii] == gc ){
      covEq = ltmp;
    }
    tmpTotal += ltmp;
  }

  (*totalCov) = tmpTotal;

  (*prop) = (double)covEq/(double)(tmpTotal);

}
/*
//discontinued
void mapOfCoverage::updateFromNoGapSamRead(pair<samRead, samRead> &pp, char *tabBases, int SHIFT){

  map<long, char> mc;
  getConsensus(mc, pp, SHIFT);

  string chrID = pp.first.rName_;

  long chrSize = this->coverage_[chrID].size() - 2;

  long pos = 0;
  char c;

  map<long, char>::iterator it;
  for(it=mc.begin() ; it!=mc.end() && pos<chrSize ; it++){
    pos = it->first;
    c = it->second;
    this->coverage_[chrID].at(pos).add(c, tabBases);
  }

}
*/

void mapOfCoverage::mergeTOgDNAmap(char *genomicMapFile, map<string, GFFtranscript> &mst, map<string, string> &mg, char *tabBases, char *tab_cpt){

  map<string, vector<char> > mStrand;
  init_mStrand(mStrand, this->coverage_);

  load_mStrand(mStrand, mst);

  mapOfCoverage mGenomic;
  mGenomic.init(mg, false);
  mGenomic.updateFromPreviousMap(genomicMapFile);

  this->mergeTOgDNAmap(mGenomic, mStrand, tabBases, tab_cpt);

}


void mapOfCoverage::mergeTOgDNAmap(mapOfCoverage &mGenomic, map<string, vector<char> > &mStrand, char *tabBases, char *tab_cpt){

  map<string, vector<coverage> >::iterator it;

  for(it=this->coverage_.begin() ; it!=this->coverage_.end() ; it++){

    string chrID = it->first;
    long lg = (it->second).size();
    for(long pos=0 ; pos<lg ; pos++){
      char strand = mStrand[chrID].at(pos);
      (it->second).at(pos).addFromGenomic( mGenomic.coverage_[chrID].at(pos), tabBases, strand, tab_cpt);
    }

  }

}




void load_mStrand(map<string, vector<char> > &mStrand, map<string, GFFtranscript> &mst){

  map<string, GFFtranscript>::iterator it;
  for(it=mst.begin() ; it!=mst.end() ; it++){
    GFFtranscript gt = it->second;
    string chrID = gt.chr_;

    char strand = gt.strand_[0];

    long start, end;
    gt.getGenomicMinMaxPos(&start, &end);

    for(long pos=start ; pos<=end ; pos++){
      char cc = mStrand[chrID].at(pos);
      if( cc == _STRAND_UNKNOWN_ ){
	mStrand[chrID].at(pos) = strand;
      } else {
	if( cc != strand ){
	  mStrand[chrID].at(pos) = _STRAND_BOTH_;
	}
      }
    }
  }

}



void init_mStrand(map<string, vector<char> > &mStrand, map<string, vector<coverage> > &mCov){

  map<string, vector<coverage> >::iterator it;
  for(it=mCov.begin() ; it!=mCov.end() ; it++){
    string chrID = it->first;
    long lg = (it->second).size();
    vector<char> v(lg, _STRAND_UNKNOWN_);
    mStrand[chrID] = v;
  }

}



/*
void mapOfCoverage::binarySave(char *fname){
  std::ofstream ofs(fname);
  boost::archive::binary_oarchive ars(ofs);
  ars & (*this);
  ofs.close();
}


void mapOfCoverage::binaryLoad(char *fname){
  std::ifstream ifs(fname);
  boost::archive::binary_iarchive arl(ifs);
  arl & (* this);
  ifs.close();
}
*/

#endif
