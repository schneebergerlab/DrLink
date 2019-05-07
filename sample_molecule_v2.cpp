/* This function aims to re-sample molecules according to base coverage distributions at each kb of given samples,
   to compare their CO frequency.

   Molecule base distributions for each kb molecule and molecule info (gzipped) from DrLink separate_kb,
   
   files will be output:

   Hequan Sun
   Date: 2018-10-07 16:30pm
*/
#include         "./gzlib/gzstream.h"
#include                     <zlib.h>
#include                   <iostream>
#include                    <sstream>
#include                    <fstream>
#include                     <vector>
#include                     <string>
#include                        <map>
#include                   <stdlib.h>
#include                    <iomanip>
#include                     <math.h>
#include                   <assert.h>
#include                     <time.h>  /* clock_t, clock, CLOCKS_PER_SEC */
#include                   <dirent.h>
#include                 <sys/stat.h>
#include                    <errno.h>
#include                   <unistd.h>
#include                  "globals.h"
#include             "split_string.h"
//
using namespace std;
//
struct MOLEINFO {
    unsigned long msta; // molecule start
    unsigned long mend; // molecule end
};
//
bool read_mcove_dist2(string                      mcfile, 
                      map<string, unsigned long>* mCov, 
                      map<string, unsigned long>* mCov_intersec,
                      map<string, int>*           mCov_intersec_samNum);             
bool write_sampled_molecule_stat2(
                      map<string, unsigned long>  molecule_cov, 
                      map<int,    unsigned long>  molecule_len, 
                      map<int,    unsigned long>  molecule_len_sum, 
                      map<int,    unsigned long>  molecule_read_num,
                      map<int,    unsigned long>  molecule_read_num_sum,
                      map<string,           int>  barcode_molecule_cnt,
                      string                      outprefix);                     
//
bool sample_molecule_v2(int argc, char* argv[])
{
    std::stringstream usage;
    usage.str("");
    if(argc < 9)
    {
        usage << endl;
        usage << "   Given sets of molecule base coverage distributions and molecule info (kb-separated), "                         
              << endl;
        usage << "   this function samples molecules according to intersection of base coverage of each kb molecule"           
              << endl;        
        const char *buildString = __DATE__", " __TIME__;
        usage << "   (compiled on " << buildString << ")"                                                   
              << endl << endl;
        usage << "   Usage: DrLink resample_v2 maxkb dataPath subfolder samples r1,r2[,...] outprefix srand-sleep" 
              << endl << endl;
        usage << "\twhere "                                                           << endl;
        usage << "\tmaxkb      : maximum molecule size to sample"                     << endl;
        usage << "\tdataPath   : path to separated molecules from DrLink separate_kb" << endl;
        usage << "\tsubfolder  : prefix to start a sample kb-molecule folder "        << endl;        
        usage << "\tsamples    : sample ids, e.g, B,C,D,E,P"                          << endl;
        usage << "\tr*         : ratio for each sample to refine base coverage "      << endl;
        usage << "\toutprefix  : general label of output files"                       << endl;
        usage << "\tsrand-sleep: integer seed value for generating random numbers"    << endl;        
        cout  << usage.str() << endl;
        return false;
    }
    double startT   = clock(); 
    // step 0. get inputs
    int maxkb = atoi(argv[2]);
    if(maxkb < 1) 
    {
        cout << "   Error: maxkb too small. Please provide a larger value w.r.t molecule size. " << endl;
        return false;
    }
    else
    {
        cout << "   Info: maximum size of molecules provided as " << maxkb << " kb. " << endl;
    }
    string dpath = (string)argv[3];
    DIR* dir = opendir(dpath.c_str());
    if (dir)
    {
        /* Directory exists. */
        closedir(dir);
    }
    else if (ENOENT == errno)
    {
        /* Directory does not exist. */
        cout << "   Error: cannot find path " << dpath << endl;
        return false;      
    }
    else ;
    string subfolder         = (string)argv[4];
    vector<string> sampid    = split_string((string)argv[5], ',');
    vector<string> readratio = split_string((string)argv[6], ',');
    assert(sampid.size()    == readratio.size());
    string outprefix         = (string)argv[7];
    int srandsleep           = atoi(argv[8]); 
    cout << "   Info: folder where kb-separate molecules are provided: "         << dpath << endl;      
    cout << "   Info: prefix of kb-molecule folder provided: "              << subfolder  << endl;
    cout << "   Info: sample ids "                                          << argv[5]    << endl;
    cout << "   Info: sample reads read ratio tuning"                       << argv[6]    << endl;                  
    cout << "   Info: sleeping " << srandsleep << " seconds to get various srand seeds. " << endl;
    
    sleep(srandsleep);
    srand (time(NULL)); 
    // prepare output files for collecting molecules & sample-wise variable for collect statistics
    string molheader = "#chr#barcode\tfirst_aligned\tlast_aligned\tmolecule_len\tmolecule_cov\tread_num\tUni_flag\tlast_aligned_end\treads_region\tR1R2\treads_id";    
    map<string, ogzstream*> allofp; 
    ////   
    map<string, map<string, unsigned long> > molecule_cov;
    map<string, map<int,    unsigned long> > molecule_len;
    map<string, map<int,    unsigned long> > molecule_len_sum;                 
    map<string, map<int,    unsigned long> > molecule_read_num;
    map<string, map<int,    unsigned long> > molecule_read_num_sum;
    map<string, map<string,           int> > barcode_molecule_cnt; 
    //
    vector<string>::iterator sitr;
    vector<string>::iterator sitr_end;
    sitr        = sampid.begin();
    sitr_end    = sampid.end();  
    int i = -1;    
    while(sitr != sitr_end)
    {
        i ++;
        if( atof(readratio[i].c_str()) > 0)
        {
            string this_id = *sitr;
            // output file
            string ofile = subfolder + "_" + this_id + "_min1000bp_sampled_molecule_table_trashme.txt.gz";
            ogzstream* f = new ogzstream(ofile.c_str(), ios::out);
            allofp.insert(std::pair<string, ogzstream*>(this_id, f) );
            cout << "   Info: output file " << ofile << " is ready. " << endl;
            *f   << molheader << endl;
            // cov 
            map<string, unsigned long> tmpcov;
            map<int,    unsigned long> tmp;        
            map<string,           int> tmpbarc;        
            molecule_cov.insert(std::pair<string, map<string, unsigned long> >(this_id, tmpcov) );
            // length
            molecule_len.insert(std::pair<string, map<int, unsigned long> >(this_id, tmp));
            // length sum
            molecule_len_sum.insert(std::pair<string, map<int, unsigned long> >(this_id, tmp));      
            // read number
            molecule_read_num.insert(std::pair<string, map<int, unsigned long> >(this_id, tmp));   
            // read number sum
            molecule_read_num_sum.insert(std::pair<string, map<int, unsigned long> >(this_id, tmp));    
            // barcode
            barcode_molecule_cnt.insert(std::pair<string, map<string, int> >(this_id, tmpbarc));
            //
        }
        sitr ++;
    }
    // step 1. sample for each kb
    for(int msize = 1; msize <= maxkb; msize ++)
    {    
        // step 1.1. get intersection value at each base coverage
        map<int, map<string, unsigned long> > mCov_individ;
        map<string, unsigned long> mCov_intersec;        
        map<string, int>           mCov_intersec_samNum; // check how many samples leading to that intersection
        int i = -1;    
        sitr        = sampid.begin();
        sitr_end    = sampid.end();            
        while(sitr != sitr_end)
        {
            string this_id = *sitr;
            map<string, unsigned long> mCov;
            std::stringstream covfile;
            covfile.str(""); 
            // dpath/reSampv8_B/B_1kb_moleCov_stat.txt
            covfile << dpath    << "/" << subfolder << "_" << this_id << "/" 
                    << this_id << "_" << msize     << "kb_moleCov_stat.txt";
            // caution: if xkb does not exist for some sample!
            if(!read_mcove_dist2(covfile.str(), &mCov, &mCov_intersec, &mCov_intersec_samNum))
            {
                cout << "   Error: reading failed on " << covfile.str() << endl;
                return false;
            }
            i ++;
            mCov_individ.insert(std::pair<int, map<string, unsigned long> >(i, mCov) );
            sitr ++;
        }
        // step 1.2. get sampling ratio
        map<int, map<string, double> > mCov_individ_sratio;
        i = -1; 
        sitr        = sampid.begin();
        sitr_end    = sampid.end();
        while(sitr != sitr_end)
        {
            string this_id = *sitr;
            map<string, double> msratio;
            i ++;
            map<string, unsigned long> mCov = mCov_individ[i];
            map<string, unsigned long>::iterator covitr;
            map<string, unsigned long>::iterator covitr_end;
            covitr     = mCov.begin();
            covitr_end = mCov.end();       
            while(covitr != covitr_end)
            {
                string key        = (*covitr).first;
                unsigned long num = (*covitr).second;      
                assert( mCov_intersec.find(key) != mCov_intersec.end() );     
                if(num == 0) num  = 1;      
                if(mCov_intersec_samNum[key]==sampid.size()) // if not all sample have this coverage, ignore it!
                {  
                    double sratio     = mCov_intersec[key]*0.80/num; // to get the minimum sampled with 0.80   
                    if(this_id.compare("B")==0 || this_id.compare("C")==0)  
                    {
                        sratio = sratio*0.90; // caution: to get a comparable molecule number to D/E/P
                    }
                    msratio.insert(std::pair<string, double>(key, sratio));
                }                                
                covitr ++;
            }
            mCov_individ_sratio.insert(std::pair<int, map<string, double> >(i, msratio));
            sitr ++;
        }
        // check
        if(true)
        {
            sitr     = sampid.begin();
            sitr_end = sampid.end();
            i = -1;   
            while(sitr != sitr_end)     
            {
                i ++;
                if( atof(readratio[i].c_str()) > 0)
                {
                    map<string, unsigned long> mCov = mCov_individ[i];
                    map<string, double>     msratio = mCov_individ_sratio[i];
                    cout << "   Check: sample " << *sitr << " msize " << msize << " kb: " << endl << endl;
                    map<string, double>::iterator sritr;
                    map<string, double>::iterator sritr_end;
                    sritr     = msratio.begin();
                    sritr_end = msratio.end();
                    while(sritr != sritr_end)
                    {
                        string key = (*sritr).first;
                        cout << "      " 
                             << (*sritr).first 
                             << "\t" << setprecision(8)  << fixed << (*sritr).second 
                             << "\tfrom\t" 
                             << mCov_intersec[key] 
                             << "/" 
                             << mCov[key] 
                             << "*0.8 (B/C-with-another0.9)\t"
                             << "sample.no=" 
                             << mCov_intersec_samNum[key]
                             << endl;
                        sritr ++;
                    }
                }
                //
                sitr ++;
            }
        }
        // step 1.3 sample molecules
        i = -1; 
        sitr        = sampid.begin();
        sitr_end    = sampid.end();
        while(sitr != sitr_end)
        {
            // get sample i
            i ++;
            map<string, double>     msratio = mCov_individ_sratio[i];
            // get sample read ratio 
            double ireadratio = atof(readratio[i].c_str());
            if(ireadratio > 0)
            {
                // get gz molecule file ready
                string this_id  = *sitr;
                ogzstream* this_f = allofp[this_id];  
                // get kb molecule file name: /dpath/reSampv8_B/B_1kb_moleCov_stat.txt
                std::stringstream kbmolfile;
                kbmolfile.str("");
                kbmolfile << dpath    << "/" << subfolder << "_" << this_id << "/" 
                        << this_id << "_" << msize     << "kb.txt";
                //
                ifstream ifp;
                ifp.open(kbmolfile.str().c_str(), ios::in);
                if(!ifp.good())
                {
                    cout << "   Error: cannot open file "       << kbmolfile.str() << endl;
                    return false;    
                }
                else
                {
                    cout    << "   Info: sampling molecules from " << kbmolfile.str() << endl;
                    //*this_f << "#   check: sampling "    << msize << " kb" << endl;
                    while(ifp.good())
                    {
                        string line("");    
                        getline(ifp, line);
                        if(line.size()==0 || line[0]=='#') continue;
                        /* get info for current molecule - 11 columns:
                           0.chr#barcode 		1#AAACACCAGAACTCGG
                           1.first_aligned 		15095089
                           2.last_aligned 		15102519
                           3.molecule_len 		7580
                           4.molecule_cov 		0.04
                           5.read_num  		 	2
                           6.Uni_flag  		 	U
                           7.last_aligned_end  		15102668
                           8.all_reads_aligned_at  	15095089,15095230;15102519,15102668
                           9.R1R2  		 	1,2
                          10.read_id 		 	rid1,rid2 
                          11.read_cigar	                cig1,cig2
                        */                    
                        vector<string> lineinfo = split_string(line, '\t');
                        if(lineinfo.size() < 12) continue;
                        double randra = rand()%1000000000*1.0 / 1000000000.0;       
                        string mckey = lineinfo[4];                    
                        if(randra <= msratio[mckey])
                        {                    
                            // stage 1. sample reads first 
                            vector<string> readPos = split_string(lineinfo[8],  ';');
                            vector<string> readR12 = split_string(lineinfo[9],  ',');
                            vector<string> readID  = split_string(lineinfo[10], ',');    
                            vector<string> readCIG = split_string(lineinfo[11], ',');
                            unsigned long first_aligned    = 0;
                            unsigned long last_aligned     = 0;
                            unsigned long mole_len         = 0;
                            unsigned long mole_covered_len = 0;        
                            double        mole_base_cov    = 0;
                            int           read_num         = 0;
                            unsigned long last_aligned_end = 0;
                            std::stringstream all_reads_aligned_at;
                            std::stringstream R1R2;
                            std::stringstream read_id;
                            std::stringstream read_cigar;
                            all_reads_aligned_at.str("");   
                            R1R2.str("");
                            read_id.str("");
                            read_cigar.str("");
                            //
                            vector<string>::iterator samread_itr;    
                            vector<string>::iterator samread_itr_end;
                            samread_itr     = readPos.begin();
                            samread_itr_end = readPos.end();    
                            int iii = -1;
                            while(samread_itr != samread_itr_end)
                            {
                                iii ++;
                                double myrand = rand()%1000000000*1.0 / 1000000000.0; 
                                //if(myrand <= readratio || iii==0 || iii==readPos.size()-1) // try keeping length of molecules
                                if(myrand <= ireadratio) // try keeping length of molecules
                                {
                                    if(false && myrand <= ireadratio)
                                    {
                                        cout << "   check: read " 
                                             << *samread_itr 
                                             << " selected, because myrand " 
                                             << myrand 
                                             << " < ireadratio " 
                                             << ireadratio << endl; 
                                    } 
                                    if(false && iii==0)
                                    {
                                        cout << "   check: read " 
                                             << *samread_itr 
                                             << " selected, because iii == 0 " 
                                             << endl;                 
                                    }    
                                    if(false && iii==readPos.size()-1)
                                    {
                                        cout << "   check: read "     
                                             << *samread_itr 
                                             << " selected, because iii ==  "  << iii << " == readPos.size()-1 == " << readPos.size()-1
                                             << endl;                 
                                    } 
                                    vector<string> readspan = split_string(*samread_itr, ',');                      
                                    if(first_aligned == 0)
                                    {
                                        first_aligned =  strtoul(readspan[0].c_str(), NULL, 0);    
                                    }
                                    last_aligned      = strtoul(readspan[0].c_str(), NULL, 0);
                                    last_aligned_end  = strtoul(readspan[1].c_str(), NULL, 0);             
                                    mole_covered_len += (last_aligned_end - last_aligned);                   
                                    //
                                    read_num ++;
                                    //
                                    if(all_reads_aligned_at.str().size() > 0)
                                    all_reads_aligned_at << ";";
                                    all_reads_aligned_at << readPos[iii];
                                    //
                                    if(R1R2.str().size()>0)
                                    R1R2                 << ",";
                                    R1R2                 << readR12[iii];
                                    //
                                    if(read_id.str().size()>0)
                                    read_id              << ",";
                                    read_id              << readID[iii];
                                    //
                                    if(read_cigar.str().size()>0)
                                    read_cigar           << ",";
                                    read_cigar           << readCIG[iii];
                                }
                                else
                                {
                                    if(false)
                                    cout << "   check: read " 
                                         << *samread_itr 
                                         << " not selected, because myrand " 
                                         << myrand 
                                         << " > ireadratio " 
                                         << ireadratio << endl;
                                }
                                samread_itr ++;
                            }        
                            mole_len      = last_aligned_end - first_aligned + 1;
                            if(mole_len < 1000) continue;
                            int    key    = (int)round((double)mole_len/1000);
                            if(key != msize) continue; // caution: molecule length might have been changed!
                            mole_base_cov = 1.0*mole_covered_len / mole_len;  
                            // collect molecule
                            *this_f<< lineinfo[0]                << "\t"
                                   << first_aligned              << "\t"
                                   << last_aligned               << "\t"    
                                   << mole_len                   << "\t"
                                   << setprecision(2)            << fixed 
                                   << mole_base_cov              << "\t"
                                   << read_num                   << "\t"
                                   << lineinfo[6]                << "\t"
                                   << last_aligned_end           << "\t"
                                   << all_reads_aligned_at.str() << "\t"
                                   << R1R2.str()                 << "\t"
                                   << read_id.str()              << "\t"
                                   << read_cigar.str()           << endl;                        
                             // collect molecule cov
                             std::stringstream updatedmckey;
                             updatedmckey.str("");
                             updatedmckey << setprecision(2) << fixed << mole_base_cov;            
                             map<string, unsigned long>::iterator mcitr = molecule_cov[this_id].find(updatedmckey.str());
                             if(mcitr == molecule_cov[this_id].end())
                             {
                                 molecule_cov[this_id].insert(std::pair<string, unsigned long>(updatedmckey.str(), 1));
                             }
                             else
                             {
                                 (*mcitr).second += 1;
                             }
                             // collect molecule len
                             map<int, unsigned long>::iterator mlitr = molecule_len[this_id].find(key);
                             if(mlitr == molecule_len[this_id].end())
                             {
                                 molecule_len[this_id].insert(std::pair<int, unsigned long>(key, 1));
                                 molecule_len_sum[this_id].insert(std::pair<int, unsigned long>(key,  mole_len));
                             }
                             else
                             {
                                 (*mlitr).second += 1;
                                 molecule_len_sum[this_id][key] += mole_len;
                             }    
                             // collect molecule read number
                             map<int,    unsigned long>::iterator readitr;
                             readitr     = molecule_read_num[this_id].find(read_num);
                             if(readitr == molecule_read_num[this_id].end())
                             {
                                 molecule_read_num[this_id].insert(std::pair<int, unsigned long>(read_num, 1));
                                 molecule_read_num_sum[this_id].insert(std::pair<int, unsigned long>(read_num, mole_covered_len));
                             }
                             else
                             {
                                 molecule_read_num[this_id][read_num] += 1;
                                 molecule_read_num_sum[this_id][read_num] += mole_covered_len;
                             }                                              
                             // molecules per barcode
                             vector<string> chrbc = split_string(lineinfo[0], '#'); // get [0]:chr, [1]:barcode/partition
                             map<string, int>::iterator bcitr;
                             bcitr = barcode_molecule_cnt[this_id].find(chrbc[1]);
                             if(bcitr == barcode_molecule_cnt[this_id].end())
                             {
                                 barcode_molecule_cnt[this_id].insert(std::pair<string, int>(chrbc[1], 1));
                             }
                             else
                             {
                                 barcode_molecule_cnt[this_id][chrbc[1]] +=1;
                             }
                        } // end of read--molecule--sampling
                    }
                }
            }
            // next sample            
            sitr ++;
        }                            
    } // end of msize
    // output statistics for each sample
    sitr        = sampid.begin();
    sitr_end    = sampid.end();
    i           = -1;
    while(sitr != sitr_end)
    {   
        i ++;
        if( atof(readratio[i].c_str()) > 0)
        { 
            string this_id  = *sitr;    
            if(!write_sampled_molecule_stat2(molecule_cov[this_id], 
                                            molecule_len[this_id], 
                                            molecule_len_sum[this_id], 
                                            molecule_read_num[this_id],
                                            molecule_read_num_sum[this_id],
                                            barcode_molecule_cnt[this_id],
                                            subfolder + "_" + this_id))
            {
                return false;
            } 
        }
        sitr ++;
    }
    // close all sampled molecule files
    map<string, ogzstream*>::iterator ofpitr;
    map<string, ogzstream*>::iterator ofpitr_end;  
    ofpitr      = allofp.begin();
    ofpitr_end  = allofp.end();    
    while(ofpitr != ofpitr_end)
    {
        ogzstream* this_f = (*ofpitr).second;
        (*this_f).close();
        ofpitr ++;
    }
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return true;
}
//
bool write_sampled_molecule_stat2(map<string, unsigned long> molecule_cov, 
                                  map<int,    unsigned long> molecule_len, 
                                  map<int,    unsigned long> molecule_len_sum, 
                                  map<int,    unsigned long> molecule_read_num,
                                  map<int,    unsigned long> molecule_read_num_sum,
                                  map<string, int>           barcode_molecule_cnt,
                                  string outprefix)
{
    // 1. molecule base coverage
    string outMCovFile = outprefix + "_min1000bp_sampled_moleCov_stat.txt";
    ofstream ofp;
    ofp.open(outMCovFile.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open " << outMCovFile << " to write data. " << endl;
        return false;
    }
    map<string, unsigned long>::iterator mcitr;
    map<string, unsigned long>::iterator mcitr_end;
    mcitr     = molecule_cov.begin();
    mcitr_end = molecule_cov.end();
    ofp << "#molecule_base_cov\tmolecule_num_in_bin" << endl;
    while(mcitr != mcitr_end)
    {
        ofp << (*mcitr).first << "\t" << (*mcitr).second << endl;
        mcitr ++;
    }
    ofp.close();
    
    // 2. molecule length: total legnth and number of molecules in bin
    string outMLenFile = outprefix + "_min1000bp_sampled_moleLen_stat.txt";    
    ofstream ofp2;
    ofp2.open(outMLenFile.c_str(), ios::out);
    if(!ofp2.good())
    {
        cout << "   Error: cannot open " << outMLenFile << " to write data. " << endl;
        return false;
    }
    map<int, unsigned long>::iterator mlitr;
    map<int, unsigned long>::iterator mlitr_end;
    mlitr     = molecule_len.begin();
    mlitr_end = molecule_len.end();
    ofp2 << "#molecule_len_kb\tmolecule_total_len_bp\tmolecule_num_in_bin" << endl;
    while(mlitr != mlitr_end)
    {
        ofp2 << (*mlitr).first  << "\t" 
             << molecule_len_sum[(*mlitr).first] << "\t" 
             << (*mlitr).second << endl;
        mlitr ++;
    }
    ofp2.close();    
    // 3. molecules per barcode/partition: find hist + output
    // find hist: 'val' (number) partitions have 'key' (number) molecules
    map<int, int> mbhist;
    map<int, int>::iterator hitr;
    map<string, int>::iterator mbcitr;
    map<string, int>::iterator mbcitr_end;
    mbcitr     = barcode_molecule_cnt.begin();
    mbcitr_end = barcode_molecule_cnt.end();
    while(mbcitr != mbcitr_end)
    {
        int cntkey = (*mbcitr).second;
        hitr = mbhist.find(cntkey);
        if(hitr == mbhist.end())
        {
            mbhist.insert(std::pair<int, int>(cntkey, 1));
        }
        else
        {
            mbhist[cntkey] += 1;
        }
        mbcitr ++;
    }  
    // output
    string outMpBHistFile = outprefix + "_min1000bp_sampled_moleNumPerBarc_stat.txt";        
    ofstream ofp3;
    ofp3.open(outMpBHistFile.c_str(), ios::out);
    if(!ofp3.good())
    {
        cout << "   Error: cannot open " << outMpBHistFile << " to write data. " << endl;
        return false;
    }   
    map<int, int>::iterator hitr_end;
    hitr     = mbhist.begin();
    hitr_end = mbhist.end();
    ofp3 << "#molecule_cnt\tpartition/barcode_num" << endl;
    while(hitr != hitr_end)
    {
        ofp3 << (*hitr).first  << "\t" << (*hitr).second << endl;
        hitr ++;
    }
    ofp3.close(); 
    // 4. number of reads per molecule
    string outReadNumFile = outprefix + "_min1000bp_sampled_readNum_stat.txt";    
    ofstream ofp4;
    ofp4.open(outReadNumFile.c_str(), ios::out);
    if(!ofp4.good())
    {
        cout << "   Error: cannot open " << outReadNumFile << " to write data. " << endl;
        return false;
    }    
    map<int,    unsigned long>::iterator readitr;
    map<int,    unsigned long>::iterator readitr_end;
    readitr     = molecule_read_num.begin();
    readitr_end = molecule_read_num.end();
    ofp4 << "#read_number\tmolecule_number_with_that_read_number\ttotal_covered_len" << endl;
    while(readitr != readitr_end)
    {
        ofp4 << (*readitr).first   << "\t" 
             << (*readitr).second  << "\t"
             << molecule_read_num_sum[(*readitr).first] << endl;
        readitr ++;
    }
    ofp4.close();
    return true;
}
//
bool read_mcove_dist2(string                     mcfile, 
                     map<string, unsigned long>* mCov, 
                     map<string, unsigned long>* mCov_intersec,
                     map<string, int>*           mCov_intersec_samNum)
{
    ifstream ifp;
    ifp.open(mcfile.c_str(), ios::in);
    if(!ifp.is_open())
    {
        cout << "   Error: cannot open file " << mcfile << endl;
        return false;
    }
    cout << "   Info: reading molecule base coverage distribution from " << mcfile << endl;
    bool first1x = true;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 2) continue;
        string        key = lineinfo[0];
        unsigned long num = strtoul(lineinfo[1].c_str(), NULL, 0);
        //
        if(atof(key.c_str()) > 1.0)
        {
            if(first1x == true)
            {
                first1x = false;
                cout << "   Warning: sample with >1.0x coverage molecules but all skipped. " << endl;
            }
            continue;
        }
        //
        map<string, unsigned long>::iterator itr;
        itr = (*mCov).find(key);
        if(itr != (*mCov).end() )
        {
            cout << "   Warning: skipping repeated line: " << line << endl;
            continue;
        }
        else
        {
            (*mCov).insert(std::pair<string, unsigned long>(key, num));
        }
        //
        itr = (*mCov_intersec).find(key);
        if(itr != (*mCov_intersec).end())
        {
            if((*itr).second > num)
            {
                // replace with lower value
                (*itr).second = num;
            }
            (*mCov_intersec_samNum)[key] += 1;
        }
        else
        {
            (*mCov_intersec).insert(std::pair<string, unsigned long>(key, num));
            (*mCov_intersec_samNum).insert(std::pair<string, int>(key, 1));
        }
    }    
    ifp.close();
    cout << "   Info: reading done. " << endl;    
    return true;
}                     
