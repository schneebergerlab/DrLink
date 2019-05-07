// globals.h
#include   <string>
#include <stddef.h>
#include      <map>
#ifndef __globals_h__
#define __globals_h__
#define INF      9999999999
#define MARGIN 0.0000000001

using namespace std;

struct BPINFO
{
    string         chr;
    unsigned long  leftp;
    unsigned long  rightp;
    unsigned long  leftMax;
    unsigned long  rightMin;
    double         score;
    string         order; // allele transition order, either "12" or "21"
    string         bpline;// raw bp info line
    string         ovline;// bp info line overlapping with the above, the smaller interval used as a representative
};

#endif
