/**

@file
@author
@date
@brief The GFFentry class represents a line of a GFF/GTF annotation file.

 */

#ifndef _GFF_HXX_
#define _GFF_HXX_


#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>



using namespace std;

class GFFentry{

public:
  /**
     Default constructor, does nothing.
   */
  GFFentry(void);
  /**
     Constructor that parses a GFF/GTF line.

     @param line A line of the GFF/GTF file
     @param sf A string containing the format of the GFF/GTF file (either GFF or GTF)
   */
  GFFentry(string line, string sf);


  string seqID_; ///< First colmn of GFF/GTF file (typically the chromosome/scaffold ID)
  string source_; ///< Second column of GFF/GTF file (the source/software used for the annotation)
  string feature_; ///< Third column of GFF/GTF file Describes the type ofentry (typically one of these: exon, intron, chromosome, scaffold, ...)

  long start_; ///< Start position of the feature (usually 1-based position from start of chromosome/scaffold)
  long end_; ///< End of the feature (usually 1-based position from start of chromosome/scaffold)
  string s_score_; ///< A score associated to the feature
  string strand_; ///< Strand on which the feature is (one of: "+", "-", or ".").

  int phase_; ///< Phase of the feature (applicable only for coding sequences/introns)

  map<string, string> mat_; ///< map of attributes

};


#endif
