#ifndef _GFFutils_HXX_
#define  _GFFutils_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>

#include "GFFentry.hxx"


#define TYPE_INTERGENIC 'T'
#define TYPE_


void getAllGFF(
	       char *fname, 
	       int nbHeaderLines, 
	       string format, 
	       map<string, vector<GFFentry> > &mvgf, 
	       vector<string>& vSources, 
	       vector<string>& vFeatures, 
	       vector<string>& vSeqNames
	       );

void getAllGFF_by_sTag(vector<GFFentry> & vgf, map<string, vector<GFFentry> > & mvgf, string sTAG, bool DEBOG=false, bool displayNBL=false);


bool my_comp_GFF(GFFentry e1, GFFentry e2);



class C_annot{
public:
  C_annot(void);
  C_annot(char strand, char type);
private:
  char strand_;
  char type_;
};



class mapOfAnnot{
public:
  mapOfAnnot(void);
  mapOfAnnot(map<string, vector<GFFentry> >);
private:
  map<string, vector<C_annot> > mAnnot_;
};


#endif
