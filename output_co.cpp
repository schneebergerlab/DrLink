#include         <map>
#include      <string>
#include      <vector>
#include   <algorithm>
#include     <iomanip>
#include     <fstream>
#include    <iostream>
#include   "globals.h"
// 
bool output_co(multimap<string, multimap<unsigned long, BPINFO> > rawbpmap, string rawbpfile)
{
    // output chr:left-coordinate-sorted raw break points: 
    // left-coordinate sorted raw break points
    fstream sortedrawfp;
    sortedrawfp.open(rawbpfile.c_str(), ios::out);
    if(!sortedrawfp.good())
    {
        cout << "   Error: cannot open file " << rawbpfile << " for collecting break point predictions. " << endl;
        return false;
    }
    sortedrawfp << "#chr\t"
                << "left_pos(sorted)\t"
                << "right_pos\t"
                << "score\t"
                << "left_alleles\t"
                << "right_alleles\t"
                << "barcode\t"
                << "vcf_mole_len\t"
                << "left_span\t"
                << "right_span\t"
                << "read_density\t"
                << "co_boundary_same_fragment\t"
                << "min_mole_RR_dist\t"
                << "max_mole_RR_dist\t"
                << "bam_mole_len\t"
                << "bam_mole_readNum\t"
                << "bam_mole_base_cov\t"          
                << "af_at_markers\t"
                << "cov_at_markers\t"
                << "co_size\t"
                << "readnum_left\t"
                << "readnum_right\t"
                << "maxium_RR_distance_in_CO\n";
                
    multimap<string, multimap<unsigned long, BPINFO> >::iterator chritr;
    multimap<string, multimap<unsigned long, BPINFO> >::iterator chritr_end;
    chritr = rawbpmap.begin();
    chritr_end = rawbpmap.end();
    while(chritr != chritr_end)
    {
        multimap<unsigned long, BPINFO> tmp;
        tmp = (*chritr).second;
        multimap<unsigned long, BPINFO>::iterator litr;
        multimap<unsigned long, BPINFO>::iterator litr_end;
        litr = tmp.begin();
        litr_end = tmp.end();
        while(litr != litr_end)
        {
            sortedrawfp << ((*litr).second).bpline << endl;
            litr ++;
        }
        chritr ++;
    }  
    sortedrawfp.close();     
    return true;
}
