/*

getIntronsFromGFF.cxx

 */


#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>


#include "GFFtranscripts.hxx"
#include "GFFtranscripts_utils.hxx"
#include "GFFutils.hxx"
#include "utils_fasta.hxx"


#define ARG_FILE_GFF 1
#define ARG_FILE_FASTA 2
#define ARG_FORMAT 3
#define ARG_sTAG 4

#define NB_BASES_IN_EXON 120
#define NB_BASES_IN_INTRON 120

using namespace std;


int main(int argc, char *argv[]){


  char tab_cpt['z'+'Z'];
  init_tab_cpt(tab_cpt);

  map<string, string> mg;
  loadFasta(mg, argv[ARG_FILE_FASTA]);

  string format(argv[ARG_FORMAT]);
  string sTAG(argv[ARG_sTAG]);


  map<string, GFFtranscript> mst;
  getAllTranscripts(argv[ARG_FILE_GFF], mst, format, sTAG, 1, mg, tab_cpt);


  fprintf(stderr, "%d transcripts ...\n", mst.size());

  map<string, GFFtranscript>::iterator it;
  for(it=mst.begin() ; it!=mst.end() ; it++){
    GFFtranscript t = it->second;
    //fprintf(stdout, ">%s (%s)\n", t.id_.c_str(), t.strand_.c_str());
    //myPrintFastaLine(t.seq_, 60, stdout);
    //fprintf(stdout, "\n");

    vector<string> vIntrons;
    t.getAllIntrons(vIntrons, mg, tab_cpt);

    vector<string> vExons;
    t.getAllExons(vExons, mg, tab_cpt);

    int nbe = vExons.size();
    for(int i=1 ; i<nbe ; i++){
      //fprintf(stdout, ">>>INTRON_%d\n", i);

      string intronSeq = vIntrons[i-1];
      string exon5p = vExons[i-1];
      string exon3p = vExons[i];

      if( intronSeq.size() > NB_BASES_IN_INTRON && exon5p.size() > NB_BASES_IN_EXON && exon3p.size() > NB_BASES_IN_EXON ){

	string b5p = exon5p.substr( exon5p.size()-NB_BASES_IN_EXON, NB_BASES_IN_EXON) + intronSeq.substr(0, NB_BASES_IN_INTRON);

	//string b3p = intronSeq.substr( intronSeq.size()-NB_BASES_IN_INTRON, NB_BASES_IN_INTRON ) + exon3p.substr(0, NB_BASES_IN_EXON);

	fprintf(
		stdout, 
		">%s_INTRON_%d_5P\n%s\n", 
		t.id_.c_str(),
		i,
		b5p.c_str()
		);
      }
    }

  }



  return 0;
}
