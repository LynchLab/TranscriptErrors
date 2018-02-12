/*

GFFutils.cxx

 */

#include <algorithm>

#include "GFFutils.hxx"


// void getAllGFF(char *fname, map<string, vector<GFFentry> > &mvgf, vector<string> vFeatures){
// Reads a GFF/GTF file and loads it into a map
// The keys of the map are the features (3rd column) in the GFF/GTF file
// To restrict the features to use: pass a vector of features in vFeatures (pass an empty vector<string> to use all features)
// Same for sources (= second column of GFF/GTF files), typically: protein_coding, ...
void getAllGFF(
	       char *fname, 
	       int nbHeaderLines, 
	       string format, 
	       map<string, vector<GFFentry> > &mvgf, 
	       vector<string> &vSources, 
	       vector<string> &vFeatures,
	       vector<string> &vSeqNames
	       ){

  string line;

  bool allFeatures = false;
  if(vFeatures.size() == 0){
    allFeatures = true;
  }

  bool allSources = false;
  if(vSources.size()==0){
    allSources = true;
  }

  bool allSeqNames = false;
  if(vSeqNames.size()==0){
    allSeqNames = true;
  }

  ifstream fic(fname);
  for(int ii=0 ; ii<nbHeaderLines ; ii++){
    getline(fic, line);
  }

  while( getline(fic, line) ){
    GFFentry gff(line, format);

    if( (allSeqNames==true || find(vSeqNames.begin(), vSeqNames.end(), gff.seqID_)!=vSeqNames.end() ) &&
	 (allFeatures==true || find(vFeatures.begin(), vFeatures.end(), gff.feature_)!=vFeatures.end() ) && 
	(allSources==true || find(vSources.begin(), vSources.end(), gff.source_)!=vSources.end() ) ){
      mvgf[gff.feature_].push_back(gff);
    }
  }

}



void getAllGFF_by_sTag(vector<GFFentry> & vgf, map<string, vector<GFFentry> > & mvgf, string sTAG, bool DEBOG, bool displayNBL){

  long nbl = vgf.size();

  for(int i=0 ; i<nbl ; i++){

    GFFentry gff = vgf[i];

    if(DEBOG==true){
      map<string, string>::iterator it;
      for(it=gff.mat_.begin() ; it!=gff.mat_.end() ; it++){
	fprintf(stdout, "MAT: %s -> %s\n", (it->first).c_str(), (it->second).c_str());
	fflush(stdout);
      }
    }

    string parentID = gff.mat_[sTAG];
    mvgf[parentID].push_back(gff);

    
    if(displayNBL==true && nbl%100000==0){
      fprintf(stderr, "%ld\n", nbl);
    }
  }

}



bool my_comp_GFF(GFFentry e1, GFFentry e2){
  if(e1.start_ < e2.start_)
    return true;
  return false;
}



