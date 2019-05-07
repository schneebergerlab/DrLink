#include         <map>
#include      <string>
#include      <vector>
#include   <algorithm>
#include     <iomanip>
#include     <fstream>
#include    <iostream>
#include   "globals.h"
//
bool postprocess_bpPrediction(multimap<string, multimap<unsigned long, BPINFO> > rawbpmap, string finalbpfile)
{
    multimap<string, multimap<unsigned long, BPINFO> > tmpmap = rawbpmap;
    // output final clusters
    fstream finfp;
    finfp.open(finalbpfile.c_str(), ios::out);
    if(!finfp.good())
    {
        cout << "   Error: cannot open file " << finalbpfile 
             << " for collecting final break point predictions. " << endl;
        return false;
    }    
    finfp << "#chr\t"
          << "left_pos\t"
          << "right_pos\t"
          << "score\t"
          << "left_alleles\t"
          << "right_alleles\t"
          << "barcode\t"
          << "molecule_len\t"
          << "left_span\t"
          << "right_span\t"
          << "Max_interval_size\t"
          << "Min_interval_size\t"
          << "support\n";
    int groupNum = 0;
    // this function overlaps predicted breakpoints
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
        // get first break point of current chr
        vector<string>        targcluster;
        vector<string>        targchr;
        vector<unsigned long> targleftp;
        vector<unsigned long> targrightp;
        vector<unsigned long> targrightMax;
        vector<unsigned long> targrightMin;
        vector<unsigned long> targleftMax;
        vector<unsigned long> targleftMin;
        vector<double>        targscore;
        vector<string>        targorder;
        targchr.push_back(((*litr).second).chr);
        targcluster.push_back(((*litr).second).bpline);
        targleftp.push_back(((*litr).second).leftp);
        targrightp.push_back(((*litr).second).rightp);
        targscore.push_back(((*litr).second).score);
        targorder.push_back(((*litr).second).order);        
        targrightMax.push_back(((*litr).second).rightp);
        targrightMin.push_back(((*litr).second).rightp);
        targleftMax.push_back(((*litr).second).leftp);
        targleftMin.push_back(((*litr).second).leftp);
        // next break point
        litr ++;
        while(litr != litr_end)
        {
            // processing here
            unsigned long ileftp  = ((*litr).second).leftp;
            unsigned long irightp = ((*litr).second).rightp; 
            string        iorder  = ((*litr).second).order;     
            
            bool grouped = false;
            for(int vi = 0; vi < targcluster.size(); vi ++)
            {
                if(iorder.compare(targorder[vi])==0 && ileftp < targrightMin[vi])
                {
                    targcluster[vi] += "\n#"+ ( (*litr).second ).bpline;
                    // right max
                    if(irightp > targrightMax[vi])
                    {
                        targrightMax[vi] = irightp;
                    }
                    // right min
                    if(irightp < targrightMin[vi])
                    {
                        targrightMin[vi] = irightp;
                    }
                    // left max
                    if(ileftp > targleftMax[vi])
                    {
                        targleftMax[vi] = ileftp;
                    }  
                    // left min
                    if(ileftp < targleftMin[vi])
                    {
                        targleftMin[vi] = ileftp;
                    }     
                    // sum of score
                    targscore[vi] += ( (*litr).second ).score;                      
                    //
                    grouped = true;
                }
            }
            if(!grouped)
            {
                targchr.push_back(((*litr).second).chr);
                targcluster.push_back(((*litr).second).bpline);
                targleftp.push_back(((*litr).second).leftp);
                targrightp.push_back(((*litr).second).rightp);
                targscore.push_back(((*litr).second).score);
                targorder.push_back(((*litr).second).order);        
                targrightMax.push_back(((*litr).second).rightp);
                targrightMin.push_back(((*litr).second).rightp);
                targleftMax.push_back(((*litr).second).leftp);
                targleftMin.push_back(((*litr).second).leftp);
            }
            
            litr ++;
        }
        // output current chr
        for(int vi = 0; vi < targcluster.size(); vi ++)
        {
            int support = std::count(targcluster[vi].begin(), targcluster[vi].end(), '#') + 1;
            // using intersection (targleftMin, targrightMax) of overlapped predictions as final interval
            finfp << "#." << ++groupNum         << "\n" 
                  << targchr[vi]                << "\t" 
                  << targleftMax[vi]            << "\t" 
                  << targrightMin[vi]           << "\t" << fixed << std::setprecision(2) 
                  << targscore[vi]/support      << "\t"
                  << targorder[vi].substr(0, 1) << "\t" 
                  << targorder[vi].substr(1, 1) << "\tBC\tmLen\tLspan\tRspan\t"
                  << targrightMax[vi] - targleftMin[vi] << "\t"
                  << targrightMin[vi] - targleftMax[vi] << "\t"
                  << support                            << "\tInterSec" << endl;
            // details of predictions      
            finfp << "#" << targcluster[vi] << endl;
            // using union (targleftMax, targrightMin) of overlapped predictions as final interval (backup)
            finfp << "##" 
                  << targchr[vi]                << "\t" 
                  << targleftMin[vi]            << "\t" 
                  << targrightMax[vi]           << "\t" << fixed << std::setprecision(2)                   
                  << targscore[vi]/support      << "\t"
                  << targorder[vi].substr(0, 1) << "\t" 
                  << targorder[vi].substr(1, 1) << "\tBC\tmLen\tLspan\tRspan\t"
                  << targrightMax[vi] - targleftMin[vi] << "\t"
                  << targrightMin[vi] - targleftMax[vi] << "\t"
                  << support                            << "\tUnion\n#" << endl;
        }
        // clear current chr
        targcluster.clear();
        targleftp.clear();
        targrightp.clear();
        targorder.clear();
        targrightMax.clear();
        targrightMin.clear();
        targleftMax.clear();
        
        chritr ++;
    }
    finfp.close();  
    return true;
}
