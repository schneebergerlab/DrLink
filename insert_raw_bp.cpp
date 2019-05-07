//
#include    <map>
#include <string>
#include "globals.h"

using namespace std;

bool insert_raw_bp(multimap<string, multimap<unsigned long, BPINFO> >* rawbpmap,
                   string        chr,
                   unsigned long leftp,
                   unsigned long rightp,
                   double        score,
                   string        order,
                   string        bpline)
{
    /*  
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
    */
    BPINFO tmpbp;
    tmpbp.chr      = chr;
    tmpbp.leftp    = leftp;
    tmpbp.rightp   = rightp;
    tmpbp.leftMax  = 0;
    tmpbp.rightMin = 0;
    tmpbp.bpline   = bpline;
    tmpbp.ovline   = "";
    tmpbp.score    = score;
    tmpbp.order    = order;
    //
    multimap<string, multimap<unsigned long, BPINFO> >::iterator chritr;
    chritr = (*rawbpmap).find(chr);
    if(chritr == (*rawbpmap).end())
    {
        // initialize map with new item of chr
        multimap<unsigned long, BPINFO> tmp;
        tmp.insert(std::pair<unsigned long, BPINFO>(leftp, tmpbp));
        //
        (*rawbpmap).insert(std::pair<string, multimap<unsigned long, BPINFO> >(chr, tmp));   
    } 
    else
    {
        // update map of chr
        ((*chritr).second).insert(std::pair<unsigned long, BPINFO>(leftp, tmpbp));
    }
    
    return true;
}
