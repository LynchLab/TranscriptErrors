/*

GFFtranscript_utils.cxx

 */


#include <algorithm>

#include "GFFtranscripts_utils.hxx"
#include "GFFutils.hxx"



void getAllTranscripts(
		       map<string, GFFtranscript> &mst, 
		       char *fname, 
		       int nbHeaderLines, 
		       string format, 
		       vector<string> vSources, 
		       string feature, 
		       string sTAG, 
		       long shift, 
		       bool addStartAndStop
		       ){

  map<string, vector<GFFentry> > mgff;

  vector<string> vFeatures;
  vFeatures.push_back(feature);
  if(addStartAndStop==true){
    vFeatures.push_back("start_codon");
    vFeatures.push_back("stop_codon");
  }

  vector<string> vSeqNames;

  getAllGFF(fname, nbHeaderLines, format, mgff, vSources, vFeatures, vSeqNames);

  map<string, vector<GFFentry> > mgff2;
  getAllGFF_by_sTag(mgff[feature], mgff2, sTAG);

  map<string, vector<GFFentry> >::iterator it;
  for(it=mgff2.begin() ; it!=mgff2.end() ; it++){
    sort( it->second.begin() , it->second.end() , my_comp_GFF );
    GFFtranscript transcript(it->second, sTAG, shift);
    mst[transcript.id_] = transcript;
  }

  if(addStartAndStop==true){
    addStartStop(mst, mgff2["stop_codon"], mgff2["stop_codon"], sTAG);
  }

}




// Simply give an empty map<string, string> to not include the sequence in the transcript
void getAllTranscripts( 
		       map<string, GFFtranscript> & mst,  
		       map<string, vector<GFFentry> > & mvGFF, 
		       string sTAG,
		       int shift,
		       map<string, string> & mg, 
		       char *tab_cpt
			){


  map<string, vector<GFFentry> >::iterator it;
  for(it=mvGFF.begin() ; it!=mvGFF.end() ; it++){
    sort( it->second.begin() , it->second.end() , my_comp_GFF );
    GFFtranscript transcript(it->second, sTAG, shift);
    mst[transcript.id_] = transcript;
  }

}



void addStartStop(map<string, GFFtranscript> &mst, vector<GFFentry> &vStarts, vector<GFFentry> &vStops, string sTag){

  vector<GFFentry>::iterator it;
  vector<GFFentry> vtab[2] = {vStarts, vStops};



  for(int i=0 ; i<2 ; i++){

    for(it=vtab[i].begin() ; it!=vtab[i].end() ; it++){

      GFFentry ge = (*it);

      long pos;

      if(ge.strand_=="+"){
	pos = ge.start_;
      } else {
	pos = ge.end_;
      }

      string parentID = ge.mat_[sTag];
      if(i==0){
	mst[parentID].posStart_ = pos;
      } else {
	mst[parentID].posStop_ = pos;
      }
    }

  }

}




void getAllTranscripts(
		       char *fname, 
		       int nbHeaderLine, 
		       map<string, GFFtranscript> &mst, 
		       string format, 
		       string sTAG, 
		       int shift, 
		       map<string, string> & mg, 
		       char *tab_cpt, 
		       vector<string> &vSources, 
		       vector<string> &vFeatures,
		       vector<string> &vSeqNames
		       ){

  // GFF entries are grouped by feature in the map mvGFFbyFeature
  map< string, vector<GFFentry> > mvGFFbyFeature;
  getAllGFF(fname, nbHeaderLine, format, mvGFFbyFeature, vSources, vFeatures, vSeqNames);


  map<string, vector<GFFentry> > mvGFF;

  map<string, vector<GFFentry> >::iterator itf;
  for(itf=mvGFFbyFeature.begin() ; itf!=mvGFFbyFeature.end() ; itf++){
    getAllGFF_by_sTag(itf->second, mvGFF, sTAG, false, false);
  }

  int nbTranscripts = mvGFF.size();

  map<string, vector<GFFentry> >::iterator itg;
  for(itg=mvGFF.begin() ; itg!=mvGFF.end() ; itg++){
    GFFtranscript gffT(itg->second, sTAG, shift, mg, tab_cpt);
    mst[gffT.id_] = gffT;
  }

}



void getAllTranscripts(char *fname, int nbHeaderLines, map<string, GFFtranscript> &mst, string format, string sTAG, int shift, map<string, string> & mg, char *tab_cpt){

  vector<string> vSources;
  vector<string> vFeatures;
  vector<string> vSeqNames;

  getAllTranscripts(fname, nbHeaderLines, mst, format, sTAG, shift, mg, tab_cpt, vSources, vFeatures, vSeqNames);

}


void getAllTranscripts(
		       char *fname, 
		       int nbHeaderLines, 
		       map<string, GFFtranscript> &mst, 
		       string format, 
		       string sTAG, 
		       int shift, 
		       vector<string> &vSources, 
		       vector<string> &vFeatures,
		       vector<string> &vSeqNames
		       ){

  map<string, string> mg;
  getAllTranscripts(fname, nbHeaderLines, mst, format, sTAG, shift, mg, NULL, vSources, vFeatures, vSeqNames);

}

void getAllTranscripts(char *fname, int nbHeaderLines, map<string, GFFtranscript> &mst, string format, string sTAG, int shift){

  vector<string> vSources;
  vector<string> vFeatures;
  vector<string> vSeqNames;

  map<string, string> mg;

  getAllTranscripts(fname, nbHeaderLines, mst, format, sTAG, shift, mg, NULL, vSources, vFeatures, vSeqNames);
}



void getPosOfIntronsGenomicScale(map<string, vector<bool> > &mpi,  map<string, GFFtranscript> &mst, map<string, string> &mg, bool exonic){

  bool intronic = true;
  if(exonic==true){ intronic = false; }

  map<string, string>::iterator it;
  for(it=mg.begin() ; it!=mg.end() ; it++){
    long lg = (it->second).size();
    vector<bool> vtmp(lg, exonic);
    mpi[(it->first)] = vtmp;
  }

  map<string, GFFtranscript>::iterator itt;

  // First, I mark as true all positions that are within an intron
  for(itt=mst.begin() ; itt!=mst.end() ; itt++){
    GFFtranscript gt = itt->second;
    gt.markIntronsGenomicScale(mpi[gt.chr_], intronic);
  }


  // Then I mark as false all the positions that are within exons (some positions can be both in intron and exon with alternative splicing, and I count those as not intronic)
  for(itt=mst.begin() ; itt!=mst.end() ; itt++){
    GFFtranscript gt = itt->second;
    gt.markExonsGenomicScale(mpi[gt.chr_], exonic);
  }
}


void transcript2genomicCoords(GFFtranscript & gt, vector<long> & vg){

  int i, e;
  int lg = 0;

  //fprintf(stderr, "Outputing exons coordinates for %s (%s):\n", gt.id_.c_str(), gt.strand_.c_str());

  int nbe = gt.vexonsCoords_.size();
  for(e=0 ; e<nbe ; e++){
    pair<long, long> pll = gt.vexonsCoords_[e];
    lg += ( (pll.second - pll.first) + 1 );
    //fprintf(stderr, "[%d - %d]\n", pll.first, pll.second);
  }
  //fprintf(stderr, "------------------------\n");


  long exonStart, exonEnd, ll;

  if( gt.strand_ == "+" ){
    for(e=0 ; e<nbe ; e++){
      exonStart = gt.vexonsCoords_[e].first;
      exonEnd = gt.vexonsCoords_[e].second;

      for(ll=exonStart ; ll<=exonEnd ; ll++){
	vg.push_back(ll);	
      }
    }
  }

  if( gt.strand_ == "-" ){
    for(e=0 ; e<nbe ; e++){
      exonStart = gt.vexonsCoords_[e].second;
      exonEnd = gt.vexonsCoords_[e].first;

      for(ll=exonStart ; ll>=exonEnd ; ll--){
	vg.push_back(ll);
      }

    }
  }


  // DEPRECATED VERSION:
  //if(gt.strand_=="+"){
  //for(e=0 ; e<nbe ; e++){
  //  pair<long, long> pll = gt.vexonsCoords_[e];
  //  for(i=pll.first ; i<=pll.second ; i++){
  //vg.push_back(i-1);
  //  }
  //}
  //}

  //if(gt.strand_=="-"){
  //for(e=nbe-1 ; e>-1 ; e--){
  //  pair<long, long> pll = gt.vexonsCoords_[e];
  //  for(i=pll.second ; i>=pll.first ; i--){
  //vg.push_back(i-1);
  //  }
  //}
  //}

}


void orderTranscriptsByChromosomeAndStartPosition(map<string, GFFtranscript> &mt, map<string, vector<GFFtranscript> > &mvt){

  map<string, GFFtranscript>::iterator itm;
  for(itm=mt.begin() ; itm!=mt.end() ; itm++){
    string chrID = (itm->second).chr_;
    mvt[chrID].push_back(itm->second);
  }

}


void getAllTranscriptsAtPosition(vector<GFFtranscript> &vtr, map<string, GFFtranscript> &mt, string & chrID, long pos){

  map<string, GFFtranscript>::iterator it;
  for(it=mt.begin() ; it!=mt.end() ; it++){
    GFFtranscript gt = it->second;
    if( gt.chr_ == chrID && gt.posMin_ <= pos && gt.posMax_ >= pos ){
      vtr.push_back(gt);
    }
  }

}


//void getAllTranscriptsAtPosition(vector<GFFtranscript> &vtr, map<string, vector<GFFtranscript> > &mvt, string chrID, long pos);

void getAllTranscriptsAtPosition(vector<GFFtranscript> &vtr, vector<GFFtranscript> &vt, long pos){

  //fprintf(stderr, "%ld\t%s\n", vt.size(), vt[0].chr_.c_str());

  vector<GFFtranscript>::iterator itv;
  for(itv=vt.begin() ; itv!=vt.end() ; itv++){
    if( itv->posMin_ <= pos && itv->posMax_>= pos ){
      vtr.push_back(*itv);
    }
  }

}


