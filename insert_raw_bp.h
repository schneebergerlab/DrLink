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
                   string        bpline);
                  
