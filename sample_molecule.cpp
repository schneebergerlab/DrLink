/* This function aims to sample molecules from pool according to intersection of their size distribution,
   to compare their CO frequency.

   Given sets of molecule size distributions and molecule info (gzipped) from DrLink molecule,
   
   files will be output:

   Hequan Sun
   Date: 2018-09-06 12:56pm
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
//<chr#bc, {msta, mend} >
multimap<string, MOLEINFO> moleMeta;
//
int srandsleep = 1;
//
bool read_msize_dist(string                      msfile, 
                     map<int, unsigned long>*    mNum, 
                     map<int, unsigned long>*    mNum_intersec);
bool read_mcove_dist(string                      mcfile, 
                     map<string, unsigned long>* mCov, 
                     map<string, unsigned long>* mCov_intersec);             
bool sampling_molecule(string                    molfile, 
                     double                      readratio,
                     string                      outprefix, 
                     map<int,           double>  msratio, 
                     map<string, unsigned long>* molecule_cov,
                     map<int,    unsigned long>* molecule_len,
                     map<int,    unsigned long>* molecule_len_sum,                      
                     map<int,    unsigned long>* molecule_read_num,
                     map<int,    unsigned long>* molecule_read_num_sum,
                     map<string,           int>* barcode_molecule_cnt);
bool write_sampled_molecule_stat(
                     map<string, unsigned long>  molecule_cov, 
                     map<int,    unsigned long>  molecule_len, 
                     map<int,    unsigned long>  molecule_len_sum, 
                     map<int,    unsigned long>  molecule_read_num,
                     map<int,    unsigned long>  molecule_read_num_sum,
                     map<string,           int>  barcode_molecule_cnt,
                     string                      outprefix);                     
//
bool sample_molecule(int argc, char* argv[])
{
    std::stringstream usage;
    usage.str("");
    if(argc < 8)
    {
        usage << endl;
        usage << "   Given sets of molecule size distributions and molecule info, "                         
              << endl;
        usage << "   this function samples molecules according to intersection of molecule sizes"           
              << endl;        
        const char *buildString = __DATE__", " __TIME__;
        usage << "   (compiled on " << buildString << ")"                                                   
              << endl << endl;
        usage << "   Usage: DrLink resample msize1.txt,msize2.txt[,...] mole1.gz,mole2.gz[,...] read1,read2[,...] labels outprefix srand-sleep-int" 
              << endl << endl;
        usage << "\twhere "                                                         << endl;
        usage << "\tmsize*.txt:   molecule size distribution of pools"              << endl;
        usage << "\tmole*.gz  :   molecule meta info corresponding to molecule size of pools from DrLink molecule"  
              << endl;
        usage << "\tread*     :   mean molecule read coverage in raw data"          << endl;
        usage << "\tlabels    :   sample-specific label of output files"            << endl;
        usage << "\toutprefix :   general label of output files"                    << endl;
        cout  << usage.str() << endl;
        return false;
    }
    double startT   = clock(); 
    int srandsleep  = atoi(argv[7]); 
    cout << "   Info: sleeping " << srandsleep << " seconds to get various srand seeds. " << endl;
    sleep(srandsleep);
    string lenfiles = (string)argv[2];
    vector<string> lenfileinfo = split_string(lenfiles, ',');
    if(lenfileinfo.size()<2)
    {
        cout << "   Error: only one file given. No need to sample!" << endl;
        return false;
    }
    string molfiles = (string)argv[3];
    vector<string> molfileinfo = split_string(molfiles, ',');
    if(molfileinfo.size()<2)
    {
        cout << "   Error: only one file given. No need to sample!" << endl;
        return false;
    }    
    if(lenfileinfo.size() != molfileinfo.size())
    {
        cout << "   Error: molecule-size file-number does not matach molecule-meta-info file-number. " << endl;
        return false;
    }
    //
    vector<string> labels = split_string((string)argv[5], ',');
    string out_prefix     = string(argv[6]);       
    // 1. read number of molecule-size distributions in pools
    map<int, unsigned long> mNum_intersec;
    map<int, map<int, unsigned long> > mNum_individ;
    vector<string>::iterator lfitr;
    vector<string>::iterator lfitr_end;
    lfitr     = lenfileinfo.begin();
    lfitr_end = lenfileinfo.end();
    int i = -1;
    while(lfitr != lfitr_end)
    {
        string msfile = *lfitr;
        map<int, unsigned long> mNum;
        if(!(read_msize_dist(msfile, &mNum, &mNum_intersec)))
        {
            cout << "   Error: failed in reading molecule size distribution. " << endl;
            return false;
        }
        i ++;
        //
        mNum_individ.insert(std::pair<int, map<int, unsigned long> >(i, mNum) );
        //
        lfitr ++;
    }
    // 2. get sampling ratio with respecitve to a specific molecule size for each pool
    map<int, map<int, double> > mNum_individ_sratio;
    lfitr     = lenfileinfo.begin();
    lfitr_end = lenfileinfo.end();
    i = -1;
    while(lfitr != lfitr_end)
    {
        map<int, double> msratio;
        //
        i ++;
        map<int, unsigned long> mNum = mNum_individ[i];
        map<int, unsigned long>::iterator lenitr;
        map<int, unsigned long>::iterator lenitr_end;
        lenitr     = mNum.begin();
        lenitr_end = mNum.end();
        while(lenitr != lenitr_end)
        {
            int key           = (*lenitr).first;
            unsigned long num = (*lenitr).second;
            assert( mNum_intersec.find(key) != mNum_intersec.end() );
            if(num == 0) num  = 1;
            double sratio     = mNum_intersec[key]*1.0/num;
            sratio = sratio*0.96;
            
            // scaling - 2018-09-25: to get C substantially sampled while keeping around 6 million molecules
            if(labels[i].compare("B") == 0) 
            {
                sratio = sratio*(6000000.0/7016664.0);
            }else
            if(labels[i].compare("C") == 0)
            {
                sratio = sratio*(6000000.0/6949772.0);
            }else
            if(labels[i].compare("D") == 0) 
            {
                sratio = sratio*(6000000.0/7273738.0);
            }else
            if(labels[i].compare("E") == 0) 
            {
                sratio = sratio*(6000000.0/7101105.0);
            }else
            if(labels[i].compare("P") == 0) 
            {
                sratio = sratio*(6000000.0/6682080.0); // 6882080
            }else ;
            //
               
            msratio.insert(std::pair<int, double>(key, sratio)); // to make D sampled
            lenitr ++;
        }
        //
        mNum_individ_sratio.insert(std::pair<int, map<int, double> >(i, msratio));      
        //
        lfitr ++; 
    }
    // check
    if(true)
    {
        lfitr     = lenfileinfo.begin();
        lfitr_end = lenfileinfo.end();
        i = -1;   
        while(lfitr != lfitr_end)     
        {
            i ++;
            map<int, unsigned long> mNum = mNum_individ[i];
            map<int, double>     msratio = mNum_individ_sratio[i];
            cout << "   Check: file " << lenfileinfo[i] << endl << endl;
            if(labels[i].compare("D") == 0) 
            {
                cout << "   Info: all length ratio decreased by 0.04 to get D a random molecule number. " << endl;
            }
            map<int, double>::iterator sritr;
            map<int, double>::iterator sritr_end;
            sritr     = msratio.begin();
            sritr_end = msratio.end();
            while(sritr != sritr_end)
            {
                int key = (*sritr).first;
                cout << "      " 
                     << (*sritr).first 
                     << "\t" << (*sritr).second 
                     << " from " 
                     << mNum_intersec[key] 
                     << "/" 
                     << mNum[key] 
                     << endl;
                sritr ++;
            }
            lfitr ++;
        }
    }
    // 3. get sampling ratio with respecitve to a specific molecule base coverage for each pool
    
    // 4. sample molecules
    // read ratio
    vector<string> readratiostr = split_string((string)argv[4], ',');
    vector<double> readratio; // ratio of reads to keep
    double minmrc = 1;
    vector<string>::iterator rritr;
    vector<string>::iterator rritr_end;
    rritr     = readratiostr.begin();
    rritr_end = readratiostr.end();
    while(rritr != rritr_end)
    {
        double rr = atof((*rritr).c_str());
        readratio.push_back(rr);
        if(minmrc > rr) minmrc = rr;
        rritr ++;
    }
    vector<double>::iterator rr2itr     = readratio.begin();
    vector<double>::iterator rr2itr_end = readratio.end();
    while(rr2itr != rr2itr_end)
    {
        (*rr2itr) = 1.0 - ((*rr2itr) - minmrc)/(*rr2itr);
        rr2itr ++;
    }
    // molecule
    vector<string>::iterator mfitr;
    vector<string>::iterator mfitr_end;
    mfitr     = molfileinfo.begin();
    mfitr_end = molfileinfo.end();
    i = -1;
    while(mfitr != mfitr_end)
    {
        i ++;
        double this_readratio       = readratio[i];
            // special - caution!
            if(labels[i].compare("B") == 0) 
            {
                this_readratio = this_readratio*0.95*0.83643;
            }else
            if(labels[i].compare("C") == 0)
            {
                this_readratio = this_readratio*1.00*0.788091;
            }else
            if(labels[i].compare("D") == 0) 
            {
                this_readratio = this_readratio*0.92*44.0/49.0;
            }else
            if(labels[i].compare("E") == 0) 
            {
                this_readratio = this_readratio*1.00;
            }else
            if(labels[i].compare("P") == 0) 
            {
                this_readratio = this_readratio*0.92*44.0/49.0;
            }else ;
        string molfile              = molfileinfo[i];
        string ilabel               = labels[i];
        map<int,    double> msratio = mNum_individ_sratio[i];    
        //    
        map<string, unsigned long> molecule_cov;
        map<int,    unsigned long> molecule_len;
        map<int,    unsigned long> molecule_len_sum;                 
        map<int,    unsigned long> molecule_read_num;
        map<int,    unsigned long> molecule_read_num_sum;
        map<string,           int> barcode_molecule_cnt;
        if(molfile.compare("skip") != 0)
        {
            if(!sampling_molecule(molfile, 
                                this_readratio,
                                out_prefix + "_" + ilabel, 
                                msratio, 
                                &molecule_cov,
                                &molecule_len,
                                &molecule_len_sum,                      
                                &molecule_read_num,
                                &molecule_read_num_sum,
                                &barcode_molecule_cnt))
            {
                return false;
            }  
            if(!write_sampled_molecule_stat(molecule_cov, 
                                            molecule_len, 
                                            molecule_len_sum, 
                                            molecule_read_num,
                                            molecule_read_num_sum,
                                            barcode_molecule_cnt,
                                            out_prefix + "_" + ilabel))
            {
                return false;
            }
        }
        mfitr ++;  
    }  
        
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return true;
}
//
bool sampling_molecule(string molfile, 
                       double readratio,
                       string outprefix, 
                       map<int,           double>  msratio, 
                       map<string, unsigned long>* molecule_cov,
                       map<int,    unsigned long>* molecule_len,
                       map<int,    unsigned long>* molecule_len_sum,                      
                       map<int,    unsigned long>* molecule_read_num,
                       map<int,    unsigned long>* molecule_read_num_sum,
                       map<string,           int>* barcode_molecule_cnt)
{
    // open for reading
    igzstream molfp;
    molfp.open(molfile.c_str());
    if(!molfp.good())
    {
        cout << "   Error: cannot file of raw molecule: " << molfile << endl;
        return false;
    }
    // open for writing
    // molecule info
    string ofilename = outprefix + "_min1000bp_sampled_molecule_table_trashme.txt.gz";
    ogzstream moutfp;
    moutfp.open(ofilename.c_str(), ios::out);
    if(!moutfp)
    {
        cout << "   Error: cannot open file " << ofilename << " for writing sampled molecules. " << endl;
        return false;
    }
    // read ids in selected molecules
    string oIDfilename = outprefix + "_min1000bp_sampled_molecule_readIDs_not_on.txt";    
    ofstream ridoutfp;
    ridoutfp.open(oIDfilename.c_str(), ios::out);
    if(!ridoutfp.is_open())
    {
        cout << "   Error: cannot open file " << oIDfilename << endl;
        return false;
    }
    //
    cout << "   Info: sampling molecules..."   << endl;
    // make sure this is out of the loop
    srand (time(NULL)); 
    while(molfp.good())
    {
        string line("");
        getline(molfp, line);
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
        //
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
            // double myrand = ((double) rand() / (RAND_MAX)); 
            double myrand = rand()%1000000000*1.0 / 1000000000.0; 
            //if(myrand <= readratio || iii==0 || iii==readPos.size()-1) // try keeping length of molecules
            if(myrand <= readratio) // try keeping length of molecules
            {
                if(false && myrand <= readratio)
                {
                    cout << "   check: read " 
                         << *samread_itr 
                         << " selected, because myrand " 
                         << myrand 
                         << " < readratio " 
                         << readratio << endl; 
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
                     << " > readratio " 
                     << readratio << endl;
            }
            samread_itr ++;
        }
        mole_len      = last_aligned_end - first_aligned + 1;
        if(mole_len < 1000) continue;
        mole_base_cov = 1.0*mole_covered_len / mole_len;
        // stage 2. sample molecule 
        
        int    key    = (int)round((double)mole_len/1000);
        // double randra = ((double) rand() / (RAND_MAX)); 
        double randra = rand()%1000000000*1.0 / 1000000000.0;         
        if(randra <= msratio[key])
        {
            // collect this molecule
            moutfp << lineinfo[0]                << "\t"
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
            std::stringstream mckey;
            mckey.str("");
            mckey << setprecision(2) << fixed << mole_base_cov;            
            map<string, unsigned long>::iterator mcitr = (*molecule_cov).find(mckey.str());
            if(mcitr == (*molecule_cov).end())
            {
                (*molecule_cov).insert(std::pair<string, unsigned long>(mckey.str(), 1));
            }
            else
            {
                (*mcitr).second += 1;
            }
            // collect molecule len
            map<int, unsigned long>::iterator mlitr = (*molecule_len).find(key);
            if(mlitr == (*molecule_len).end())
            {
                (*molecule_len).insert(std::pair<int, unsigned long>(key, 1));
                (*molecule_len_sum).insert(std::pair<int, unsigned long>(key,  mole_len));
            }
            else
            {
                (*mlitr).second += 1;
                (*molecule_len_sum)[key] += mole_len;
            }
            // collect molecule read number
            map<int,    unsigned long>::iterator readitr;
            readitr     = (*molecule_read_num).find(read_num);
            if(readitr == (*molecule_read_num).end())
            {
                (*molecule_read_num).insert(std::pair<int, unsigned long>(read_num, 1));
                (*molecule_read_num_sum).insert(std::pair<int, unsigned long>(read_num, mole_covered_len));
            }
            else
            {
                (*molecule_read_num)[read_num] += 1;
                (*molecule_read_num_sum)[read_num] += mole_covered_len;
            }
            // molecules per barcode
            vector<string> chrbc = split_string(lineinfo[0], '#'); // get [0]:chr, [1]:barcode/partition
            map<string, int>::iterator bcitr;
            bcitr = (*barcode_molecule_cnt).find(chrbc[1]);
            if(bcitr == (*barcode_molecule_cnt).end())
            {
                (*barcode_molecule_cnt).insert(std::pair<string, int>(chrbc[1], 1));
            }
            else
            {
                (*barcode_molecule_cnt)[chrbc[1]] +=1;
            }     
            // read ids
            if(false)
            {
                vector<string> readids = split_string(read_id.str(), ',');
                map<string, int> checkreads;
                vector<string>::iterator iditr;
                vector<string>::iterator iditr_end;
                iditr     = readids.begin();
                iditr_end = readids.end();
                while(iditr != iditr_end)
                {
                    string this_id = *iditr;
                    if(checkreads.find(this_id) == checkreads.end() )
                    {
                        checkreads.insert(std::pair<string, int>(this_id, 1) );
                        ridoutfp << this_id << endl;
                    }
                    iditr ++;
                }
                checkreads.clear();
            }
        }
    }
    molfp.close();  
    moutfp.close();    
    ridoutfp.close();  
    return true;
}
//
bool write_sampled_molecule_stat(map<string, unsigned long> molecule_cov, 
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
bool read_msize_dist(string msfile, map<int, unsigned long>* mNum, map<int, unsigned long>* mNum_intersec)
{
    ifstream ifp;
    ifp.open(msfile.c_str(), ios::in);
    if(!ifp.is_open())
    {
        cout << "   Error: cannot open file " << msfile << endl;
        return false;
    }
    cout << "   Info: reading molecule size distribution from " << msfile << endl;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 2) continue;
        int           key = atoi(lineinfo[0].c_str());
        unsigned long len = strtoul(lineinfo[lineinfo.size()-1].c_str(), NULL, 0);
        //
        map<int, unsigned long>::iterator itr;
        itr = (*mNum).find(key);
        if(itr != (*mNum).end() )
        {
            cout << "   Warning: skipping repeated line: " << line << endl;
            continue;
        }
        else
        {
            (*mNum).insert(std::pair<int, unsigned long>(key, len));
        }
        //
        itr = (*mNum_intersec).find(key);
        if(itr != (*mNum_intersec).end())
        {
            if((*itr).second > len)
            {
                // replace with lower value
                (*itr).second = len;
            }
        }
        else
        {
            (*mNum_intersec).insert(std::pair<int, unsigned long>(key, len));
        }
    }
    ifp.close();
    cout << "   Info: reading done. " << endl;
    return true;
}
//
bool read_mcove_dist(string                      mcfile, 
                     map<string, unsigned long>* mCov, 
                     map<string, unsigned long>* mCov_intersec)
{
    ifstream ifp;
    ifp.open(mcfile.c_str(), ios::in);
    if(!ifp.is_open())
    {
        cout << "   Error: cannot open file " << mcfile << endl;
        return false;
    }
    cout << "   Info: reading molecule base coverage distribution from " << mcfile << endl;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 2) continue;
        std::stringstream key0;
        key0.str(""); 
        key0 << lineinfo[0];
        string        key = key0.str();
        unsigned long num = strtoul(lineinfo[1].c_str(), NULL, 0);
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
        }
        else
        {
            (*mCov_intersec).insert(std::pair<string, unsigned long>(key, num));
        }
    }    
    ifp.close();
    cout << "   Info: reading done. " << endl;    
    return true;
}                     
