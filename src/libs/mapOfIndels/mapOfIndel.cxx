/*

mapOfIndel.cxx

 */


#include "remappedRead/remappedRead_utils.hxx"
#include "TRerrors/TRerrors.hxx"

#include "mapOfIndel.hxx"

//#define MIN_POS_READ 2
//#define MAX_POS_READ 103
//#define MIN_POS_MATE 0
//#define MAX_POS_MATE 99


using namespace std;

C_mapOfIndel::C_mapOfIndel(void){
  ;
}



C_mapOfIndel::C_mapOfIndel(map<string, string> &mg, bool verbose){

  this->init(mg, verbose);

}


void C_mapOfIndel::init(map<string, string> &mg, bool verbose){

  map<string, string>::iterator it;

  for(it=mg.begin() ; it!=mg.end() ; it++){
    string chrID = it->first;
    long lg = (it->second).size();

    if( verbose==true ){
      fprintf(stderr, ", %s", chrID.c_str());
    }

    vector<C_indelCoverage> vtmp(lg);
    this->indelCoverage_[chrID] = vtmp;
  }

}

/*
void C_mapOfIndel::updateFromRemapped(char *fname, map<string, string> &mBarcodes, bool verbose){
 
 map<string, vector< pair<remappedRead, remappedRead> > > m;
  bool includeSequence = false;
  bool computeSeqWithIndels = false;
  bool includeQuality = false;
  initMapGroupedPairedReadsBy_mRNA(fname, m, '0', 'z', includeSequence, computeSeqWithIndels, includeQuality);
  this->updateFromRemapped(m, 
			   mBarcodes, 
			   verbose);
  m.clear();

}


void C_mapOfIndel::updateFromRemapped(
				    map<string, vector< pair<remappedRead, remappedRead> > > & mr,
				    map<string, string> &mBarcodes,
				    bool verbose
				    ){


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

	this->updateFromFamily(itb->second);

      }

    }

  }

}
*/


void C_mapOfIndel::updateFromFamily( vector< pair<remappedRead, remappedRead> > &vFamily ){

  int i;
  string chrID = vFamily[0].first.rName_;

  map<long, char> mConsensusIndel;
  getMapConsensusIndel(mConsensusIndel, vFamily);

  map<long, char>::iterator it;
  for(it=mConsensusIndel.begin() ; it!=mConsensusIndel.end() ; it++){
    long pos = it->first;    
    char stat = it->second;
    if(stat!='U'){
      this->indelCoverage_[chrID].at(pos).cov_++;
    }

    if(stat=='I'){
      this->indelCoverage_[chrID].at(pos).ins_++;
    }

    if(stat=='D'){
      this->indelCoverage_[chrID].at(pos).del_++;
    }

  }

}




void C_mapOfIndel::updateFromPreviousMap(char *fname){
  
  string line, chrID, stmp;
  long pos, cov, ins, del;

  ifstream fic(fname);
  while( getline(fic, line) ){
    istringstream iss(line);
    getline(iss, chrID, '\t');

    getline(iss, stmp, '\t'); // position
    pos = atol(stmp.c_str());

    getline(iss, stmp, '\t');
    cov = atol(stmp.c_str());

    getline(iss, stmp, '\t');
    ins = atol(stmp.c_str());

    getline(iss, stmp, '\t');
    del = atol(stmp.c_str());

    this->indelCoverage_[chrID].at(pos).cov_ += cov;
    this->indelCoverage_[chrID].at(pos).ins_ += ins;
    this->indelCoverage_[chrID].at(pos).del_ += del;

  }

}



void C_mapOfIndel::display(FILE *nf){

  map<string, vector<C_indelCoverage> >::iterator it;

  for(it=this->indelCoverage_.begin() ; it!=this->indelCoverage_.end() ; it++){
    string chrID = it->first;
    long lg = (it->second).size();
    for(long pos=0 ; pos<lg ; pos++){
      fprintf(nf, "%s\t%d\t", chrID.c_str(), pos);
      (it->second).at(pos).display(nf);
      fprintf(nf, "\n");
    }
  }

}


int getListOfIndels(
		     vector< pair<remappedRead, remappedRead> > &vr, 
		     vector< pair<long, long> > &vIns,
		     vector< pair<long, long> > &vDel,
		     bool overlapOnly
		     ){

  if(vr.size()<1){
    return 0;
  }

  int nbIns = vr[0].first.vInsertions_.size();
  int nbDel = vr[0].first.vDeletions_.size();

  for(int ii=0 ; ii<nbIns ; ii++){
    pair<long, long> ppi = vr[0].first.vInsertions_[ii];
    if( vr[0].second.asInsertion(ppi)==true || (vr[0].second.coversPosition(ppi.first)==false && overlapOnly==false) ){
      vIns.push_back(ppi);
    }
  }

  for(int id=0 ; id<nbDel ; id++){
    pair<long, long> ppd = vr[0].first.vDeletions_[id];
    if( vr[0].second.asDeletion(ppd)==true || (vr[0].second.coversPosition(ppd.first)==false && overlapOnly==false) ){
      vDel.push_back(ppd);
    }
  }

  nbIns = vr[0].second.vInsertions_.size();
  nbDel = vr[0].second.vDeletions_.size();

  for(int ii=0 ; ii<nbIns ; ii++){
    pair<long, long> ppi = vr[0].second.vInsertions_[ii];
    if( vr[0].first.coversPosition(ppi.first)==false && overlapOnly==false ){
      vIns.push_back(ppi);
    }
  }

  for(int id=0 ; id<nbDel ; id++){
    pair<long, long> ppd = vr[0].second.vDeletions_[id];
    if( vr[0].first.coversPosition(ppd.first)==false &&overlapOnly==false ){
      vDel.push_back(ppd);
    }
  }

  int nbIndel = vIns.size() + vDel.size();
  if( nbIndel>1 ){ return nbIndel; }

  // Now I make sure that the indel is present in all members of the round:
  int sizeRound = vr.size();

  for(int rr=1 ; rr<sizeRound ; rr++){

    vector< pair<long, long> >::iterator iti;
    for(iti=vIns.begin() ; iti!=vIns.end() ; iti++){
      pair<long, long> ppi = (*iti);
      if( (vr[rr].first.coversPosition(ppi.first)==true && vr[rr].first.asInsertion(ppi)==false) ||
	  (vr[rr].second.coversPosition(ppi.first)==true && vr[rr].second.asInsertion(ppi)==false)
	  ){
	vIns.erase(iti);
	iti--;
      }
    }

    vector< pair<long, long> >::iterator itd;
    for(itd=vDel.begin() ; itd!=vDel.end() ; itd++){
      pair<long, long> ppd = (*itd);
      if( ( vr[rr].first.coversPosition(ppd.first)==true && vr[rr].first.asDeletion(ppd)==false ) ||
	  ( vr[rr].second.coversPosition(ppd.first)==true && vr[rr].second.asDeletion(ppd)==false )
	  ){
	vDel.erase(itd);
	itd--;
      }
    }
  }

  nbIns = vIns.size();
  nbDel = vDel.size();
  nbIndel = nbIns + nbDel;

  return nbIndel;
}

void getMapConsensusIndel(map<long, char> &m, vector< pair<remappedRead, remappedRead> > &vFamily){

  vector< pair<long, long> > vIns;
  vector< pair<long, long> > vDel;

  bool onlyOverlap = false;
  int nbIndel = getListOfIndels(vFamily, vIns, vDel, onlyOverlap);
  int nbIns = vIns.size();
  int nbDel = vDel.size();

  if( nbIndel > 1 ){
    return;
  }

  int sizeMap = vFamily[0].first.vmaps_.size();
  for(int ii=0 ; ii<sizeMap ; ii++){
    long start = vFamily[0].first.vmaps_[ii].first;
    long end = (start + vFamily[0].first.vmaps_[ii].second) - 1;
    for(long pos=start ; pos<=end ; pos++){
      m[pos] = 'C';
    }
  }

  sizeMap = vFamily[0].second.vmaps_.size();
  for(int ii=0 ; ii<sizeMap ; ii++){
    long start = vFamily[0].second.vmaps_[ii].first;
    long end = (start + vFamily[0].second.vmaps_[ii].second) - 1;
    for(long pos=start ; pos<=end ; pos++){
      m[pos] = 'C';
    }
  }

  for(int ii=0 ; ii<nbIns ; ii++){
    pair<long, long> ppi = vIns[ii];
    long start = ppi.first;
    long end = (start + ppi.second) - 1;
    for(long pos=start ; pos<=end ; pos++){
      m[pos] = 'I';
    }
  }

  for(int id=0 ; id<nbDel ; id++){
    pair<long, long> ppd = vDel[id];
    long start = ppd.first;
    long end = (start + ppd.second) - 1;
    for(long pos=start ; pos<=end ; pos++){
      m[pos] = 'D';
    }
  }

}




int updateIndelErrors(
		       C_mapOfIndel &moi,
		       int minCov,
		       double maxProp,
		       vector< vector< pair<remappedRead, remappedRead> > > &vFam,
		       FILE *nfCandidates,
		       long *NB_INS,
		       long *NB_DEL,
		       long *NB_POS_COVERED_FOR_INDELS,
		       map<string, vector<int> > &mhp,
		       map<int, pair<long, long> > &mObsHP,
		       map<string, GFFtranscript> &mGFF,
		       map<string, vector<long> > &mObs,
		       map<string, string> &mg,
		       char *tab_cpt,
		       char *tabBases
		      ){

  
  int i;
  int NB_ROUND = vFam.size();


  int MIN_POS_READ = 2;
  int MAX_POS_READ = 103;
  int MIN_POS_MATE = 0;
  int MAX_POS_MATE = 99;

  int runSize = -1;
 
  vector< pair<long, long> > vIns;
  vector< pair<long, long> > vDel;
  pair<remappedRead, remappedRead> prr;
  getSharedIndels(vFam, NB_ROUND, vIns, vDel, prr);


  int readMapSize = prr.first.vmaps_.size();
  int mateMapSize = prr.second.vmaps_.size();

  int overSeq;

  if( prr.first.strand_ == "+" ){
    overSeq = prr.first.vmaps_[0].first - prr.second.vmaps_[0].first;
  } else {
    overSeq = prr.second.vmaps_[mateMapSize-1].first - prr.first.vmaps_[readMapSize-1].first;
  }

  if( overSeq > 0 ){
    MIN_POS_MATE = overSeq + 4;
    //fprintf(stderr, "Mate goes past read (overSeq=%d) : %s\n", overSeq, prr.first.qName_.c_str());
    //return 0;
  }

  /*
  if( (prr.first.strand_=="+" && prr.second.vmaps_[0].first < prr.first.vmaps_[0].first) || 
      (prr.first.strand_=="-" && prr.second.vmaps_[mateMapSize-1].first > prr.first.vmaps_[readMapSize-1].first )
      ){
    fprintf(stderr, "Mate goes past read : %s\n", prr.first.qName_.c_str());
    return 0;
  }
  */

  
  int nbIns = vIns.size();
  int nbDel = vDel.size();
  
  //fprintf(stderr, "INS:%d\tDEL:%d\n", nbIns, nbDel);
  
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // I refuse all families with more than one indel!
  if( (nbIns + nbDel) > 1 ){
    return 0;
  }
  
  char gc;

  string chrID = prr.first.rName_;
  string tID = prr.first.transcriptID_;
  bool cpt = false;
  if( prr.first.strand_=="-" ){ cpt = true; }

  GFFtranscript gff = mGFF[tID];
  
  //map<long, pair<char, char> > mIntersection;
  //getIntersection(prr.first, prr.second, mIntersection, tab_cpt);
  
  map<long, char> mhc;
  map<long, char> msc;
  getConsensus(vFam, mhc, msc, tabBases, tab_cpt, false);

  if( msc.size() == 0 ){
    //fprintf(stderr, "msc.size()==0 for %s\n", prr.first.qName_.c_str());
    //getConsensus(vFam, mhc, msc, tabBases, tab_cpt, false);
    return 0;
  }


  int nbMismatchesInSoftAln = 0;
  int posTmp = 0;
  map<long, char>::iterator it;

  for(it=msc.begin() ,posTmp=0 ; it!=msc.end() ; it++, posTmp++){
    long gPos = it->first;
    gc = mg[chrID][gPos];
    if(cpt==true){ gc = tab_cpt[gc]; }
    if( gc != it->second ){
      nbMismatchesInSoftAln++;
    }
  }

  if( nbMismatchesInSoftAln > 1 ){
    return 0;
  }

  //fprintf(stderr, "MSC:%d:%d\n", msc.size(), nbMismatchesInSoftAln);

  // Going through the covered positions
  for(it=msc.begin() ; it!=msc.end() ; it++){

    long gPos = it->first;

    long readCoord = prr.first.genomeToReadCoords(gPos);
    long mateCoord = prr.second.genomeToReadCoords(gPos);

    if( (readCoord<MIN_POS_READ || readCoord>MAX_POS_READ) && (mateCoord<MIN_POS_MATE || mateCoord>MAX_POS_MATE) )
      continue;


    C_indelCoverage covi = moi.indelCoverage_[chrID].at(gPos);
    long cov = covi.cov_;
    long ins = covi.ins_;
    long del = covi.del_;
    long indel = ins + del;
    double prop = (double)indel/(double)cov;

    if(minCov==0 || maxProp==0.0 || (cov>=minCov && prop<=maxProp) ){
      mObs[chrID].at(gPos)++;
      (*NB_POS_COVERED_FOR_INDELS)++;
      long tPos = gff.genomicToSplicedTranscript(gPos);

      if( tPos>0 && tPos<mhp[tID].size() ){
	runSize = mhp[tID].at(tPos);
	(mObsHP[runSize].second)++;
      } else {
	runSize = -1;
	fprintf(stderr, "WARNING: tPos=%d (for gPos=%d and tID=%s)\n", tPos, gPos, tID.c_str());
      }
    }
  }


  // Calling all the insertions:
  for(i=0 ; i<nbIns ; i++){
    long gPos = vIns[i].first;
    long indelSize = vIns[i].second;

    long readCoord = prr.first.genomeToReadCoords(gPos);
    long mateCoord = prr.second.genomeToReadCoords(gPos);

    if( (readCoord<MIN_POS_READ || readCoord>MAX_POS_READ) && (mateCoord<MIN_POS_MATE || mateCoord>MAX_POS_MATE) ){
      return 0;
    }

    remappedRead tmpRead;
    long coord;
    if( readCoord>0 ){ 
      tmpRead = prr.first;
      coord = readCoord;
    } else {
      tmpRead = prr.second;
      coord = mateCoord;
    }
    //tmpRead.getSeqWithIndels();
    //fprintf(stderr, "Extracting at position (%d , %d) from %s\n", coord, indelSize, tmpRead.seqWithIndels_.c_str());
    //fflush(stderr);
    string indelSeq = tmpRead.seq_.substr(coord, indelSize);


    C_indelCoverage covi = moi.indelCoverage_[chrID].at(gPos);
    long cov = covi.cov_;
    long ins = covi.ins_;
    long del = covi.del_;
    long indel = ins + del;
    double prop = (double)indel/(double)cov;

    if(minCov==0 || maxProp==0.0 || (cov>=minCov && prop<=maxProp) ){

      long tPos = gff.genomicToSplicedTranscript(gPos);
      if( tPos>0 && tPos<mhp[tID].size() ){
	runSize = mhp[tID].at(tPos);
	(mObsHP[runSize].first)++;
      } else {
	runSize = -1;
	fprintf(stderr, "WARNING: tPos=%d (for gPos=%d and tID=%s)\n", tPos, gPos, tID.c_str());
      }

      (*NB_INS)++;
      //fprintf(stderr, "strand : %s\n", prr.first.strand_.c_str());
      displayIndelCandidate(nfCandidates, vFam, 'I', indelSize, indelSeq, chrID, gPos, tID, tPos, runSize, prr.first.strand_, cov, ins, del, prop);
      //fprintf(nfCandidates, "msc.size() = %d\n", msc.size());
      //displayMSC(nfCandidates, msc, mg, chrID, cpt, tab_cpt);
      fflush(nfCandidates);
    }
  }


  // Calling the deletions:
  for(i=0 ; i<nbDel ; i++){
    long gPos = vDel[i].first;
    long indelSize = vDel[i].second;

    long readCoord = prr.first.genomeToReadCoords(gPos);
    long mateCoord = prr.second.genomeToReadCoords(gPos);

    if( (readCoord<MIN_POS_READ || readCoord>MAX_POS_READ) && (mateCoord<MIN_POS_MATE || mateCoord>MAX_POS_MATE) ){
      return 0;
    }

    string indelSeq("-", indelSize);

    C_indelCoverage covi = moi.indelCoverage_[chrID].at(gPos);
    long cov = covi.cov_;
    long ins = covi.ins_;
    long del = covi.del_;
    long indel = ins + del;
    double prop = (double)indel/(double)cov;
    
    if(minCov==0 || maxProp==0.0 || (cov>=minCov && prop<=maxProp) ){

      long tPos = gff.genomicToSplicedTranscript(gPos);
      if( tPos>0 && tPos<mhp[tID].size() ){
	runSize = mhp[tID].at(tPos);
	(mObsHP[runSize].first)++;
      } else {
	runSize = -1;
	fprintf(stderr, "WARNING: tPos=%d (for gPos=%d and tID=%s)\n", tPos, gPos, tID.c_str());
      }

      (*NB_DEL)++;
      //fprintf(stderr, "strand : %s\n", prr.first.strand_.c_str());
      displayIndelCandidate(nfCandidates, vFam, 'D', indelSize, indelSeq, chrID, gPos, tID, tPos, runSize, prr.first.strand_, cov, ins, del, prop);
      //fprintf(nfCandidates, "msc.size() = %d\n", msc.size());
      //displayMSC(nfCandidates, msc, mg, chrID, cpt, tab_cpt);
      //fflush(nfCandidates);
    }    
  }

  return 1;

}


void displayIndelCandidate(FILE *nf, vector< vector< pair<remappedRead, remappedRead> > >&vFam, char stat, int indelSize, string indel, string chrID, long gPos, string tID, long tPos, int runSize, string strand, long cov, long ins, long del, double prop){

  fprintf(nf, 
	  "%c\t%d\t%s\t",
	  stat,
	  indelSize,
	  indel.c_str()
	  );


  fprintf(nf,
	  "%s\t%d\t%s\t%d\t%d\t%s\t",
	  chrID.c_str(),
	  gPos,
	  tID.c_str(),
	  tPos,
	  runSize,
	  strand.c_str()
	  );

  fprintf(nf,
	  "cov:%d\tins:%d\tdel:%d\tprop:%.4f",
	  cov,
	  ins,
	  del,
	  prop
	  );


  int famSize = vFam.size();
  for(int i=0 ; i<famSize ; i++){
    fprintf(nf, "\t{");
    for(int j=0 ; j<vFam[i].size() ; j++){
      if(j>0){
	fprintf(nf, " , ");
      }

      fprintf(nf, "%s", vFam[i].at(j).first.qName_.c_str());
    }
    fprintf(nf, "}");
  }

  fprintf(nf, "\n");

}


void getSharedIndels(vector< vector< pair<remappedRead, remappedRead> > > &vFam, 
		     int NB_ROUND,
		     vector< pair<long, long> > &vIns, 
		     vector< pair<long, long> > &vDel,
		     pair<remappedRead, remappedRead> &prrArg
		     ){

  int iRound, roundSize;

  vector< pair<long, long> > vInsTmp;
  vector< pair<long, long> > vDelTmp;

  pair<remappedRead, remappedRead> pp;
  bool ppInitialized = false;

  for(iRound=0 ; iRound<NB_ROUND ; iRound++){
    roundSize = vFam[iRound].size();
    for(int jj=0 ; jj<roundSize ; jj++){
      if( ppInitialized == false){
	pair<remappedRead, remappedRead> prrTmp = vFam[iRound].at(jj);
	prrArg = prrTmp;
	initIndelsTmp(prrTmp, vInsTmp, vDelTmp );
	ppInitialized = true;
      } else {
	updateIndelsTmp( vFam[iRound].at(jj), vInsTmp, vDelTmp );
      }
    }
  }

  int nbIns = vInsTmp.size();
  int nbDel = vDelTmp.size();

  for(int i=0 ; i<nbIns ; i++){
    vIns.push_back(vInsTmp[i]);
  }

  for(int i=0 ; i<nbDel ; i++){
    vDel.push_back(vDelTmp[i]);
  }

}


void initIndelsTmp(pair<remappedRead, remappedRead> &prr, vector< pair<long, long> > &vIns, vector< pair<long, long> > &vDel){

  int i;
  int nbIns = prr.first.vInsertions_.size();
  for(i=0 ; i<nbIns ; i++){
    pair<long, long> ppIns = prr.first.vInsertions_[i];
    if( prr.second.coversPosition(ppIns.first)==false || prr.second.asInsertion(ppIns)==true ){
      vIns.push_back(ppIns);
    }
  }
  nbIns = prr.second.vInsertions_.size();
  for(i=0 ; i<nbIns ; i++){
    pair<long, long> ppIns = prr.second.vInsertions_[i];
    if( prr.first.coversPosition(ppIns.first)==false ){
      vIns.push_back(ppIns);
    }
  }


  int nbDel = prr.first.vDeletions_.size();
  for(i=0 ; i<nbDel ; i++){
    pair<long, long> ppDel = prr.first.vDeletions_[i];
    if( prr.second.coversPosition(ppDel.first)==false || prr.second.asDeletion(ppDel)==true ){
      vDel.push_back(ppDel);
    }
  }
  nbDel = prr.second.vDeletions_.size();
  for(i=0 ; i<nbDel ; i++){
    pair<long, long> ppDel = prr.second.vDeletions_[i];
    if( prr.first.coversPosition(ppDel.first)==false ){
      vDel.push_back(ppDel);
    }
  }

}


void updateIndelsTmp(pair<remappedRead, remappedRead> &prr, vector< pair<long, long> > &vIns, vector< pair<long, long> > &vDel){

  vector< pair<long, long> >::iterator it;
  for(it=vIns.begin() ; it!=vIns.end() ; it++){
    if( ( prr.first.coversPosition(it->first)==true && prr.first.asInsertion((*it))==false ) || 
	( prr.second.coversPosition(it->first)==true && prr.second.asInsertion((*it))==false )
	){
      vIns.erase(it);
      it--;
    }
  }

  for(it=vDel.begin() ; it!=vDel.end() ; it++){
    if( ( prr.first.coversPosition(it->first)==true && prr.first.asDeletion((*it))==false ) || 
	( prr.second.coversPosition(it->first)==true && prr.second.asDeletion((*it))==false )
	){
      vDel.erase(it);
      it--;
    }
  }

}


void displayMSC(FILE *nf, map<long, char> &msc, map<string, string> &mg, string &chrID, bool cpt, char *tab_cpt){

  char gc;
  map<long, char>::iterator it;

  for(it=msc.begin() ; it!=msc.end() ; it++){
    long gPos = it->first;
    gc = mg[chrID][gPos];
    if(cpt==true){ gc = tab_cpt[gc]; }

    fprintf(nf, "%d\t%c:%c\n", gPos, gc, it->second);
    
  }

}



int C_mapOfIndel::updateFromFamily(
				   vector< vector< pair<remappedRead, remappedRead> > > &vFam,
				   FILE *nfCandidates,
				   long *NB_INS,
				   long *NB_DEL,
				   long *NB_POS_COVERED_FOR_INDELS,
				   map<string, string> &mg,
				   char *tab_cpt,
				   char *tabBases
				   ){

  
  int i;
  int NB_ROUND = vFam.size();

  int MIN_POS_READ = 2;
  int MAX_POS_READ = 103;
  int MIN_POS_MATE = 0;
  int MAX_POS_MATE = 99;

  int runSize = -1;
  long tPos = -1;
 
  vector< pair<long, long> > vIns;
  vector< pair<long, long> > vDel;
  pair<remappedRead, remappedRead> prr;
  getSharedIndels(vFam, NB_ROUND, vIns, vDel, prr);


  int readMapSize = prr.first.vmaps_.size();
  int mateMapSize = prr.second.vmaps_.size();

  int overSeq;

  if( prr.first.strand_ == "+" ){
    overSeq = prr.first.vmaps_[0].first - prr.second.vmaps_[0].first;
  } else {
    overSeq = prr.second.vmaps_[mateMapSize-1].first - prr.first.vmaps_[readMapSize-1].first;
  }

  if( overSeq > 0 ){
    MIN_POS_MATE = overSeq + 4;
    //fprintf(stderr, "Mate goes past read (overSeq=%d) : %s\n", overSeq, prr.first.qName_.c_str());
    //return 0;
  }

  int nbIns = vIns.size();
  int nbDel = vDel.size();
  
  //fprintf(stderr, "INS:%d\tDEL:%d\n", nbIns, nbDel);
  
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // I refuse all families with more than one indel!
  if( (nbIns + nbDel) > 1 ){
    return 0;
  }
  
  char gc;
  string chrID = prr.first.rName_;
  string tID = prr.first.transcriptID_;
  bool cpt = false;
  if( prr.first.strand_=="-" ){ cpt = true; }
  
  map<long, char> mhc;
  map<long, char> msc;
  getConsensus(vFam, mhc, msc, tabBases, tab_cpt, false);

  if( msc.size() == 0 ){
    //fprintf(stderr, "msc.size()==0 for %s\n", prr.first.qName_.c_str());
    //getConsensus(vFam, mhc, msc, tabBases, tab_cpt, false);
    return 0;
  }


  int nbMismatchesInSoftAln = 0;
  int posTmp = 0;
  map<long, char>::iterator it;

  for(it=msc.begin() ,posTmp=0 ; it!=msc.end() ; it++, posTmp++){
    long gPos = it->first;
    gc = mg[chrID][gPos];
    if(cpt==true){ gc = tab_cpt[gc]; }
    if( gc != it->second ){
      nbMismatchesInSoftAln++;
    }
  }

  if( nbMismatchesInSoftAln > 2 ){
    return 0;
  }

  // Going through the covered positions
  for(it=msc.begin() ; it!=msc.end() ; it++){

    long gPos = it->first;

    long readCoord = prr.first.genomeToReadCoords(gPos);
    long mateCoord = prr.second.genomeToReadCoords(gPos);

    if( (readCoord<MIN_POS_READ || readCoord>MAX_POS_READ) && (mateCoord<MIN_POS_MATE || mateCoord>MAX_POS_MATE) )
      continue;

    this->indelCoverage_[chrID].at(gPos).cov_++;

    (*NB_POS_COVERED_FOR_INDELS)++;

  }


  // Calling all the insertions:
  for(i=0 ; i<nbIns ; i++){
    long gPos = vIns[i].first;
    long indelSize = vIns[i].second;
    
    long readCoord = prr.first.genomeToReadCoords(gPos);
    long mateCoord = prr.second.genomeToReadCoords(gPos);
    
    if( (readCoord<MIN_POS_READ || readCoord>MAX_POS_READ) && (mateCoord<MIN_POS_MATE || mateCoord>MAX_POS_MATE) ){
      return 0;
    }
    
    int incre = 1;
    if(cpt==true){
      incre = -1;
    }

    long gPosTmp;
    for(int ii=0, gPosTmp=gPos ; ii<indelSize ; ii++, gPosTmp+=incre){
      this->indelCoverage_[chrID].at(gPosTmp).ins_++;
    }

    (*NB_INS)++;
    //fprintf(stderr, "strand : %s\n", prr.first.strand_.c_str());
    displayIndelCandidate(nfCandidates, vFam, 'I', indelSize, chrID, gPos, tID, tPos, runSize, prr.first.strand_, readCoord, mateCoord);
    //fprintf(nfCandidates, "msc.size() = %d\n", msc.size());
    //displayMSC(nfCandidates, msc, mg, chrID, cpt, tab_cpt);
    fflush(nfCandidates);

  }


  // Calling the deletions:
  for(i=0 ; i<nbDel ; i++){
    long gPos = vDel[i].first;
    long indelSize = vDel[i].second;

    long readCoord = prr.first.genomeToReadCoords(gPos);
    long mateCoord = prr.second.genomeToReadCoords(gPos);

    if( (readCoord<MIN_POS_READ || readCoord>MAX_POS_READ) && (mateCoord<MIN_POS_MATE || mateCoord>MAX_POS_MATE) ){
      return 0;
    }


    int incre = 1;
    if(cpt==true){
      incre = -1;
    }

    long gPosTmp;
    for(int ii=0, gPosTmp=gPos ; ii<indelSize ; ii++, gPosTmp+=incre){
      this->indelCoverage_[chrID].at(gPosTmp).del_++;
    }


    (*NB_DEL)++;
    //fprintf(stderr, "strand : %s\n", prr.first.strand_.c_str());
    displayIndelCandidate(nfCandidates, vFam, 'D', indelSize, chrID, gPos, tID, tPos, runSize, prr.first.strand_, readCoord, mateCoord);
    //fprintf(nfCandidates, "msc.size() = %d\n", msc.size());
    //displayMSC(nfCandidates, msc, mg, chrID, cpt, tab_cpt);
    //fflush(nfCandidates);
    
  }

  return 1;
  
}



void displayIndelCandidate(FILE *nf, vector< vector< pair<remappedRead, remappedRead> > >&vFam, char stat, int indelSize, string chrID, long gPos, string tID, long tPos, int runSize, string strand, long readCoord, long mateCoord){

  string indel = "";
  long cov = 0;
  long ins = 0;
  long del = 0;
  double prop = 0.0;

  displayIndelCandidate(nf, vFam, stat, indelSize, indel, chrID, gPos, tID, tPos, runSize, strand, cov, ins, del, prop);
  fprintf(nf, "readCoord:%d\tmateCoord:%d\n", readCoord, mateCoord);

}
