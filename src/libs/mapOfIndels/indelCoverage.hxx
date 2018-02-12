#ifndef _indelCoverage_hxx_
#define _indelCoverage_hxx_


#include <stdio.h>

#include <vector>
#include <map>
#include <iostream>

using namespace std;

class C_indelCoverage{

public:
  C_indelCoverage(void);
  C_indelCoverage(long cov, long ins, long del);

  void update(long cov, long ins, long del);
  void display(FILE *nf);

  long cov_;
  long ins_;
  long del_;
};


#endif
