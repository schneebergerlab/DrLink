/* This function samples molecules to compare relative CO frequency in different samples.

   Given out_pool_*_25000_min1000bp_molecule_table_trashme.txt.gz of several samples,
   
   for each sample:
	=>1. build molecule size distribution: 1kb:n1, 2kb:n2, ..., 65kb:n65
	=>2. for each size, build base coverage disitriubtion: 0.01:a1, 0.02:a2, ..., 20:a200
	   

   Hequan Sun
   Date: 2018-10-06 14:19pm
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
bool write_molecule_kbcov(map<string, unsigned long> molecule_kbcov, 
                     string molfolder, 
                     string mollabel, int mlkey);
//
bool separate_kb_molecule(int argc, char* argv[])
{
    std::stringstream usage;
    usage.str("");
    if(argc < 5)
    {
        usage << endl;
        usage << "   Given sets of molecule info (from DrLink molecule), "                         
              << endl;
        usage << "   this function samples molecules according to intersection of molecule sizes and base coverage"           
              << endl;
        const char *buildString = __DATE__", " __TIME__;
        usage << "   (compiled on " << buildString << ")"                                                   
              << endl << endl;
        usage << "   Usage: DrLink separate_kb mole1.gz,mole2.gz[,...] label1,label2[,...] outprefix" 
              << endl << endl;
        usage << "\twhere "                                                         << endl;
        usage << "\tmole*.gz  :   molecule meta info corresponding to molecule size of pools from DrLink molecule"  
              << endl;
        usage << "\tlabel*    :   sample-specific label of output files"            << endl;
        usage << "\toutprefix :   general label of output files"                    << endl;
        cout  << usage.str() << endl;
        return false;
    }
    double startT   = clock(); 
    // get inputs
    vector<string> mfileinfo  = split_string((string)argv[2], ',');
    vector<string> mlabelinfo = split_string((string)argv[3], ',');
    if(mfileinfo.size() != mlabelinfo.size())
    {
        cout << "   Error: number of molecule files does not match number of labels. Exited." << endl;
        return false;
    }
    string outprefix = (string)argv[4];
    if(outprefix.size()==0)
    {
        outprefix = "Enjo";
    }
    //
    vector<string>::iterator mfitr;
    vector<string>::iterator mfitr_end;
    mfitr     = mfileinfo.begin();
    mfitr_end = mfileinfo.end();
    int i = -1;
    string molheader = "#chr#barcode\tfirst_aligned\tlast_aligned\tmolecule_len\tmolecule_cov\tread_num\tUni_flag\tlast_aligned_end\treads_region\tR1R2\treads_id";
    map<int, map<string,  unsigned long> > molecule_cov; // map<kb, map<0.01x,0.02x,..., number of molecules> >
    while(mfitr != mfitr_end)
    {
        i ++;
        string molfile  = *mfitr;        
        string mollabel = mlabelinfo[i];
        cout << "   Info: building statistics for " << mollabel << ":" << molfile << endl;
        unsigned long raw_mole_num      = 0;
        unsigned long selected_mole_num = 0;
        // create an intermediate folder for collecting details about a CO-molecule
        string molfolder = outprefix+"_"+mollabel;
        DIR* dir = opendir(molfolder.c_str());
        if (dir)
        {
            /* Directory exists. */
            closedir(dir);
        }
        else if (ENOENT == errno)
        {
            /* Directory does not exist. */
            const int dir_err = mkdir(molfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (dir_err == -1)
            {
                cout << "   Error: cannot create directory " << molfolder << endl;
                return false;
            }        
        }
        else ;
        // read raw molecule meta info file
        igzstream molfp;
        molfp.open(molfile.c_str());
        if(!molfp.good())
        {
            cout << "   Error: cannot open file of raw molecule: " << molfile << endl;
            return false;
        }
        while(molfp.good())
        {
            string line("");
            getline(molfp, line);
            if(line.size()==0 || line[0]=='#') continue;            
            vector<string> lineinfo = split_string(line, '\t');
            raw_mole_num ++;
            // molecule size
            int mlen = atoi(lineinfo[3].c_str());
            if(mlen < 1000) continue;
            int mlkey = (int)round((double)mlen/1000);
            map<int, map<string, unsigned long> >::iterator mlitr;
            mlitr = molecule_cov.find(mlkey); // to check if key exists
            if(mlitr == molecule_cov.end())
            {
                map<string, unsigned long> tmpcovmap;
                molecule_cov.insert(std::pair<int, map<string, unsigned long> >(mlkey, tmpcovmap) );
            }
            mlitr = molecule_cov.find(mlkey); // to update key: value
            // molecule cov
            string mckey = lineinfo[4];
            map<string, unsigned long>::iterator mcitr = ((*mlitr).second).find(mckey);
            if(mcitr == (*mlitr).second.end())
            {
                ((*mlitr).second).insert(std::pair<string, unsigned long>(mckey, 1));
            }
            else
            {
                (*mcitr).second += 1;
            }
            // output x.kb molecule to x.kb file
            std::stringstream kbmfile;
            kbmfile.str(""); 
            kbmfile << "./" << molfolder << "/" << mollabel << "_" << mlkey << "kb.txt";
            ofstream outkbfp;
            outkbfp.open(kbmfile.str().c_str(), ios::out | ios::app);
            if(!outkbfp.good())
            {
                cout   << "   Error: cannot open file " << kbmfile.str() << endl;
                return false;
            }
            else
            {
                outkbfp << line << endl;
                outkbfp.close();
                selected_mole_num ++;
            }
            //
            if(raw_mole_num%1000000 == 0)
            {
                cout << "   Info: " 
                     << raw_mole_num 
                     << " raw molecules, where " 
                     << selected_mole_num 
                     << " selected. " 
                     << endl;
            }
        }
        // output sample kb molecule base coverage
        map<int, map<string, unsigned long> >::iterator kbitr;
        map<int, map<string, unsigned long> >::iterator kbitr_end;
        kbitr     = molecule_cov.begin();
        kbitr_end = molecule_cov.end();
        while(kbitr != kbitr_end)
        {
            if(!write_molecule_kbcov((*kbitr).second, molfolder, mollabel, (*kbitr).first))
            {
                cout << "   Error: cannot write molecule cov for " << (*kbitr).first << " kb." << endl;
                return false;
            }
            kbitr ++;
        }
        // next sample
        mfitr ++;
    }
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return true;
}
//
bool write_molecule_kbcov(map<string, unsigned long> molecule_kbcov, string molfolder, string mollabel, int mlkey)
{
    // 1. molecule base coverage
    std::stringstream outMCovFile;
    outMCovFile.str(""); 
    outMCovFile << "./" << molfolder << "/" << mollabel << "_" << mlkey << "kb_moleCov_stat.txt";
    ofstream ofp;
    ofp.open( (outMCovFile.str()).c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open " << outMCovFile.str() << " to write data. " << endl;
        return false;
    }
    map<string, unsigned long>::iterator mcitr;
    map<string, unsigned long>::iterator mcitr_end;
    mcitr     = molecule_kbcov.begin();
    mcitr_end = molecule_kbcov.end();
    ofp << "#molecule_base_cov\tmolecule_num_in_bin" << endl;
    while(mcitr != mcitr_end)
    {
        ofp << (*mcitr).first << "\t" << (*mcitr).second << endl;
        mcitr ++;
    }
    ofp.close();
    
    return true;
}
