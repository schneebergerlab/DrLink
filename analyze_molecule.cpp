/* this function analyze/cluster barcodes for regions of a chr to get

   1) distribution of molecule length              - x: molecule_len                      , y: molecule_num
   2) distribution of molecule base coverage       - x: read_num * read_len / molecule_len, y: molecule_num
   3) number of moleculer per partition (optional) - x: molecule_num                      , y: partition_num
   4) number of reads per molecule                 - x: read_num                          , y: molecule_num
   
   input: a 8 column gzipped file converted from bam in a previous step.

*/
#include   "./gzlib/gzstream.h"
#include            <iostream>
#include             <sstream>
#include             <fstream>
#include              <vector>
#include              <string>
#include                 <map>
#include            <stdlib.h>
#include             <iomanip>
#include              <math.h>
#include            <assert.h>
#include              <time.h>  /* clock_t, clock, CLOCKS_PER_SEC */

#include "split_string.h"
//#include "read_chrsize.h"

using namespace std;
//
struct ALIGNINFO {
    unsigned long first;  // alignment start
    unsigned long second; // alignment end
    string        rorder; // read ordering in pairs: R1 or R2
    string        readid; // read id
    string        rcigar; // CIGAR string for read
};
// barcode statistics: TODO
struct BCINFO {
    int           molecnt; // total number of molecules
    int           molelen; // total length of molecules
    int           readcnt; // total number of reads
    unsigned long basecnt; // total number of bases
};
//
int  decipher_cigar(string cigar);   
bool write_molecule_stat(map<string, unsigned long> molecule_cov, 
                         map<int,    unsigned long> molecule_len, 
                         map<int,    unsigned long> molecule_len_sum, 
                         map<int,    unsigned long> molecule_read_num,
                         map<int,    unsigned long> molecule_read_num_sum,
                         map<int,    unsigned long> molecule_read_interval,
                         map<string, int>           barcode_molecule_cnt,
                         string out_prefix,
                         int minMoleSize);
bool write_chr_alignment(ogzstream& moutfp,
                         multimap<string, ALIGNINFO>* bcPosQueue, 
                         map<string,  unsigned long>* molecule_cov,
                         map<int,     unsigned long>* molecule_len,
                         map<int,     unsigned long>* molecule_len_sum, 
                         map<int,     unsigned long>* molecule_read_num,
                         map<int,     unsigned long>* molecule_read_num_sum,
                         map<int,     unsigned long>* molecule_read_interval,
                         map<string,  int>*           barcode_molecule_cnt,
                         string lastChr,
                         unsigned long winsize,
                         int minMoleSize);                                                   

bool analyze_molecule(int argc, char* argv[])
{
    if(argc != 6)
    {
        cout << "\n   Usage: DrLink molecule preprocessed_bamInfo.gz max_inter_read_dist min_molecule_size out_prefix_str "               << endl;
        cout << "     max_inter_read_dist\tINT\tbp: reads with same barcode but distance >= value will be separated into sub-molecules. " << endl;
        cout << "                             \tThis depends on sequencing library design: if expected mean size of molecules is EXP kp, we can use 1.5*EXP bp." << endl;   
        cout << "     min_molecule_size\tINT\tbp: size of molecules <= value will not contribute to final statistics, e.g., 1000 bp. "    << endl;
        return false;
    }
    clock_t tbeg;
    tbeg = clock();
    bool verbose = true;            // TODO: add option
    unsigned long winsize = 100000; // same barcode in chr with distance>100kb => two molecules; note this depends on sequencing library design
    // input file
    string procbamfile("");    
    int minMoleSize = 1000;
    std::stringstream minMsize;
    minMsize.str("");
    minMsize << minMoleSize;    
    string out_prefix("");
    if(argc == 6)
    {
        procbamfile = (string)argv[2];
        winsize     = strtoul(argv[3], NULL, 0); 
        minMoleSize = atoi(argv[4]);
        out_prefix  = (string)argv[5]; 
        minMsize.str("");
        minMsize << minMoleSize;       
    }
    // step 1. read chr size info
    /*
    map<string, string> chrSizes;
    if(!read_chrsize(chrsizefile.c_str(), &chrSizes))
    {
        cout << "   Error: reading chr size failed (analyze_molecule.cpp). " << endl;
        return false;
    }
    */
    // step 2.
    cout << "   Info: analyzing molecules started..." << endl;
    // open converted file.gz.
    if(procbamfile.find(".gz")==std::string::npos)
    {
        cout << "   Error: current version only supported gzipped file. " << endl;
        return false;
    }
    igzstream alifp;
    alifp.open(procbamfile.c_str());
    if(!alifp.good())
    {
        cout << "   Error: cannot open " << procbamfile << endl;
        return false;
    }
    // 2.1 prepare variable for analyzing/outputing alignments/clusters
    string ofilename = out_prefix + "_min" + minMsize.str() + "bp_molecule_table_trashme.txt.gz";
    ogzstream moutfp;
    moutfp.open(ofilename.c_str(), ios::out);
    if(!moutfp)
    {
        cout << "    Error: cannot open file " << ofilename << " for writing. " << endl;
        return false;
    }
    moutfp << "#chr#barcode\tfirst_aligned\tlast_aligned\tmolecule_len\tmolecule_cov\tread_num\tUni_flag\tlast_aligned_end\treads_region\tR1R2\treads_id" << endl;
    map<string,  unsigned long> molecule_cov;          // key: 0.01x, 0.02x,0.03x, ...,    val: number of molecules
    map<int,     unsigned long> molecule_len;          // key: 0kb, 1kb, 2kb,...           val: number of molecules.
    map<int,     unsigned long> molecule_len_sum;      // key: 0kb, 1kb, 2kb,...           val: total length of molecules in bin 
    map<int,     unsigned long> molecule_read_num;     // key: 1,2,3...                    val: number of molecules
    map<int,     unsigned long> molecule_read_num_sum; // key: 1,2,3...                    val: total length of reads in bin
    // note that below is a statistic of the read-distances within the same chr#barcode (but possibly several molecules).
    map<int,     unsigned long> molecule_read_interval;// key: 200,201,...                 val: total number of cases with that inter-read-distance-given-by-key
    map<string,            int> barcode_molecule_cnt;  // key: barcode1, barcode2,...      val: number of molecules
    multimap<string, ALIGNINFO> bcPosQueue; // from barcode find read-aligning (sta, end)-positions
    // 2.2 cluster existing alignments in queue and output molecule info -- caution in RAM
    string lastChr("");
    string currentChr("");
    unsigned long lineNum = 0;
    while(alifp.good())
    {
        string line("");
        getline(alifp, line);
        if(line.size()==0 || line[0]=='#' || line[0]=='*') continue;
        
        lineNum ++;
        if(lineNum%10000000 == 0) cout << "    " << lineNum << "th aligm..." << endl;
        
        // chr pos inbamCIGAR barcode MI read-ordering:inbamFLAG reversed-align read-id
        // 1	125	128M	TTTCCTCAGGGTTTCT	0	R1:99 rc        idddddd
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 6) continue;        
        if(lastChr.size() == 0)
        {
            lastChr    = lineinfo[0];
        }
        currentChr = lineinfo[0];        
        
        if(currentChr.compare(lastChr)!=0 && bcPosQueue.size()>0)
        {
            // output alignments for last chr and record statistics
            write_chr_alignment(moutfp,
                                &bcPosQueue, 
                                &molecule_cov,
                                &molecule_len,
                                &molecule_len_sum,
                                &molecule_read_num,
                                &molecule_read_num_sum,
                                &molecule_read_interval,
                                &barcode_molecule_cnt,
                                lastChr,
                                winsize,
                                minMoleSize);
            // clean info for last Chr
            bcPosQueue.clear();
            lastChr = currentChr;
        }
        unsigned long stapos = strtoul(lineinfo[1].c_str(), NULL, 0); // alignment start position
        int covlen = 0;
        covlen = decipher_cigar(lineinfo[2]);                         // covered len by alignment
        string key("");
        key    = lineinfo[0] + "#" + lineinfo[3];                     // chr+"#"+barcode            
        ALIGNINFO span;
        span.first  = stapos;
        span.second = stapos+covlen-1;                                //
        if(lineinfo[5].find("R1") != std::string::npos)
        {
            span.rorder = "1";
        }
        else
        if(lineinfo[5].find("R2") != std::string::npos)
        {
            span.rorder = "2";
        }  
        else
        {
            span.rorder = "x";
        }      
        span.readid = lineinfo[7];
        span.rcigar = lineinfo[2];
        // chr+"#"+barcode => (sta, end) : allow multi-aligns with the same barcode.
        bcPosQueue.insert(std::pair<string, ALIGNINFO>(key, span));
    }
    // output alignments for last chr record statistics
    if(bcPosQueue.size()>0)
    {
        // output alignments for last chr and record statistics
        write_chr_alignment(moutfp,
                            &bcPosQueue, 
                            &molecule_cov,
                            &molecule_len,
                            &molecule_len_sum,
                            &molecule_read_num,       
                            &molecule_read_num_sum, 
                            &molecule_read_interval,                    
                            &barcode_molecule_cnt,
                            lastChr,
                            winsize,
                            minMoleSize);
        // clean info for last Chr
        bcPosQueue.clear();
    }
    //
    alifp.close();
    moutfp.close();    
    // output statistics
    if(!write_molecule_stat(molecule_cov, 
                            molecule_len, 
                            molecule_len_sum,
                            molecule_read_num,
                            molecule_read_num_sum,
                            molecule_read_interval,
                            barcode_molecule_cnt, 
                            out_prefix, 
                            minMoleSize))
    {
        return false;
    }
    //
    cout << "   Info: molecule raw info    collected in " << ofilename << endl;
    cout << "   Info: molecule coverage    collected in"  << out_prefix + "_min" + minMsize.str() + "bp_moleCov_stat.txt" << endl;
    cout << "   Info: molecule length      collected in"  << out_prefix + "_min" + minMsize.str() + "bp_moleLen_stat.txt" << endl;  
    cout << "   Info: molecule read number collected in"  << out_prefix + "_min" + minMsize.str() + "bp_readNum_stat.txt" << endl;  
    cout << "   Info: reads' distances     collected in " << out_prefix + "_min" + minMsize.str() + "bp_readDist_stat.txt" << endl;
    cout << "   Info: time on analyzing molecules: "      << (float)(clock()-tbeg)/CLOCKS_PER_SEC << " seconds.\n" << endl;        
    return true;
}

bool write_chr_alignment(ogzstream& moutfp,
                         multimap<string, ALIGNINFO>* bcPosQueue, 
                         map<string, unsigned long>* molecule_cov,
                         map<int,    unsigned long>* molecule_len,
                         map<int,    unsigned long>* molecule_len_sum, 
                         map<int,    unsigned long>* molecule_read_num,
                         map<int,    unsigned long>* molecule_read_num_sum,
                         map<int,    unsigned long>* molecule_read_interval,
                         map<string, int>*           barcode_molecule_cnt,
                         string lastChr,
                         unsigned long winsize,
                         int minMoleSize)
{
    // bcPosQueue: multimap<chr#barcode, <read.align.start, read.align.end> >
    multimap<string, ALIGNINFO>::iterator mitr; 
    multimap<string, ALIGNINFO>::iterator mitr_end;
    mitr     = (*bcPosQueue).begin();
    mitr_end = (*bcPosQueue).end();
    long moleNum  = 0;
    long alignNum = 0;
    while(mitr != mitr_end)
    {
        string key = (*mitr).first;
        vector<string> chrbc = split_string(key, '#'); // get [0]:chr, [1]:barcode/partition
        std::pair <multimap<string, ALIGNINFO>::iterator, 
                   multimap<string, ALIGNINFO>::iterator> clitr_region;
        clitr_region = (*bcPosQueue).equal_range(key);
        multimap<string, ALIGNINFO>::iterator clitr;
        unsigned long firstAlignSta  = 0; // for calculating molecule length; molecule first alignment       
        unsigned long firstAlignEnd  = 0; // not really used
        unsigned long lastAlignSta   = 0; // molecule last alignment       
        unsigned long lastAlignEnd   = 0; // for calculating molecule length; for breaking molecule
        unsigned long lastR1AlignSta = 0; // for calculating R1 distance
        unsigned long lastR1AlignEnd = 0; // ..
        unsigned long lastR2AlignSta = 0; // for calculating R2 distance
        unsigned long lastR2AlignEnd = 0; // ..     
        int     readNum = 0;
        int     covleng = 0;
        int moleculeLen = 0;
        std::stringstream alignedPos;     // aligned positions of reads of a molecule
        alignedPos.str("");
        //
        std::stringstream readids; // id of reads
        readids.str("");
        std::stringstream rcigars; // cigar of reads
        rcigars.str("");
        std::stringstream rpairids;// id of R1/R2
        rpairids.str("");
        //
        map<int, unsigned long> im_read_interval;
        for (clitr=clitr_region.first; clitr!=clitr_region.second; ++clitr) // span of reads within chr#barcode
        {
            ALIGNINFO span         = (*clitr).second;       
            if(clitr==clitr_region.first || (span.first>lastAlignEnd && span.first-lastAlignEnd+1>winsize)) // new molecule
            {
                // summary last cluster with the same barcode
                if(firstAlignSta != 0)
                {
                   assert(lastAlignEnd>firstAlignSta);
                   moleculeLen = lastAlignEnd-firstAlignSta+1; // in bp
                   if(moleculeLen==0) moleculeLen = 1; // caution: avoid potential dividing 0
                   moutfp << lastAlignSta     << "\t" 
                          << moleculeLen      << "\t" 
                          << setprecision(2)  << fixed << (float)covleng/moleculeLen << "\t" 
                          << readNum          << "\tR\t" 
                          << lastAlignEnd << "\t"
                          << alignedPos.str() << "\t" 
                          << rpairids.str()   << "\t" 
                          << readids.str()    << "\t"
                          << rcigars.str()    << endl;                          
                   if(moleculeLen >= minMoleSize)
                   {
                       // molecule cov
                       std::stringstream mckey;
                       mckey.str("");
                       mckey << setprecision(2) << fixed << (float)covleng/moleculeLen;
                       map<string, unsigned long>::iterator mcitr = (*molecule_cov).find(mckey.str());
                       if(mcitr == (*molecule_cov).end())
                       {
                           (*molecule_cov).insert(std::pair<string, unsigned long>(mckey.str(), 1));
                       }
                       else
                       {
                           (*mcitr).second += 1;
                       }
                       // molecule len
                       int mlkey = (int)round((double)moleculeLen/1000); 
                       map<int, unsigned long>::iterator mlitr = (*molecule_len).find(mlkey);
                       if(mlitr == (*molecule_len).end())
                       {
                           (*molecule_len).insert(std::pair<int, unsigned long>(mlkey, 1));
                           (*molecule_len_sum).insert(std::pair<int, unsigned long>(mlkey, moleculeLen));
                       }
                       else
                       {
                           (*mlitr).second += 1;
                           (*molecule_len_sum)[mlkey] += moleculeLen;
                       }
                       // molecules per barcode
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
                       // molecule read number
                       map<int,    unsigned long>::iterator readitr;
                       readitr = (*molecule_read_num).find(readNum);
                       if(readitr == (*molecule_read_num).end())
                       {
                           (*molecule_read_num).insert(std::pair<int, unsigned long>(readNum, 1));
                           (*molecule_read_num_sum).insert(std::pair<int, unsigned long>(readNum, span.second - span.first + 1));
                       }
                       else
                       {
                           (*molecule_read_num)[readNum] += 1;
                           (*molecule_read_num_sum)[readNum] += span.second - span.first + 1;
                       }
                       // molecule read distances:
                       map<int, unsigned long>::iterator imrditr;
                       map<int, unsigned long>::iterator imrditr_end;
                       map<int, unsigned long>::iterator imrditr_final;
                       imrditr     = im_read_interval.begin();
                       imrditr_end = im_read_interval.end();
                       while(imrditr != imrditr_end)
                       {
                           //if(abs((*imrditr).first) > 500) // note: with 0.1x of 50kb molecule, dist between 18 readpairs is 2777bp
                           if(true)
                           {
                               imrditr_final = (*molecule_read_interval).find((*imrditr).first);
                               if(imrditr_final == (*molecule_read_interval).end())
                               {
                                   (*molecule_read_interval).insert(std::pair<int, unsigned long>((*imrditr).first, (*imrditr).second));
                               }
                               else
                               {
                                   (*imrditr_final).second += (*imrditr).second;
                               }
                           }
                           imrditr ++;
                       }        
                   }
                   //
                   readNum = 0;
                   // aligned positions of reads of a molecule
                   alignedPos.str("");
                   // read id and order of pair info
                   readids.str("");
                   rcigars.str("");
                   rpairids.str("");                   
                   covleng = 0;
                }
                // start new/next cluster
                moleNum ++;
                firstAlignSta = span.first;
                firstAlignEnd = span.second;
                moutfp << (*clitr).first << "\t" << firstAlignSta << "\t";
                // read dist info: if do following, then distance of reads between inter-molecules is ignored.
                lastR1AlignSta = 0; // initializing
                lastR1AlignEnd = 0;  
                lastR2AlignSta = 0; // initializing
                lastR2AlignEnd = 0;    
                im_read_interval.clear();             
            }
            // get covered length and read number
            covleng += span.second - span.first + 1;        
            readNum ++;
            // aligned positions of reads of a molecule: each item is a spanning range of a read
            if(alignedPos.str().size()>0)
            {
                alignedPos  << ";" << span.first << "," << span.second;
            }
            else
            {
                alignedPos  << span.first << "," << span.second;
            }
            // read id and order of pair info
            if(readids.str().size()>0)
            {
                readids  << "," << span.readid;
                rcigars  << "," << span.rcigar;
                rpairids << "," << span.rorder; 
            }
            else
            {
                readids  << span.readid;
                rcigars  << span.rcigar;
                rpairids << span.rorder; 
            }            
            alignNum ++;
            // get R1.vs.R1-/R2.vs.R2-read distances within the chr#barcode
            if(span.rorder.compare("1")==0)
            {
                if(lastR1AlignEnd != 0)// if there is one R1-read before current R1-read
                {
                    map<int, unsigned long>::iterator mritr;
                    int distSize = 0;
                    if(span.first>=lastR1AlignEnd)
                    {
                        distSize = span.first-lastR1AlignEnd+1;
                    }
                    else
                    {
                        distSize = -1*(lastR1AlignEnd - span.first+1); // in case overlapping
                    }
                    mritr = im_read_interval.find(distSize);
                    if(mritr == im_read_interval.end())
                    {
                        im_read_interval.insert(std::pair<int, unsigned long>(distSize, 1));
                    }
                    else
                    {
                        (*mritr).second += 1;
                    }
                    //cout << "   check: " << (*clitr).first << ": span.first-lastR1AlignEnd=abs(" << span.first << "-" << lastR1AlignEnd << ")=" << distSize << endl;
                }
                lastR1AlignSta = span.first; // updating
                lastR1AlignEnd = span.second;
            }
            else
            if(span.rorder.compare("2")==0)
            {
                if(lastR2AlignEnd != 0) // if there is one R2-read before current R2-read
                {
                    map<int, unsigned long>::iterator mritr;
                    int distSize = 0;
                    if(span.first>=lastR2AlignEnd)
                    {
                        distSize = span.first-lastR2AlignEnd+1;
                    }
                    else
                    {
                        distSize = -1*(lastR2AlignEnd - span.first+1); // in case overlapping
                    }
                    mritr = im_read_interval.find(distSize);
                    if(mritr == im_read_interval.end())
                    {
                        im_read_interval.insert(std::pair<int, unsigned long>(distSize, 1));
                    }
                    else
                    {
                        (*mritr).second += 1;
                    }
                    //cout << "   check: " << (*clitr).first << ": span.first-lastR2AlignEnd=abs(" << span.first << "-" << lastR2AlignEnd << ")=" << distSize << endl;
                }
                lastR2AlignSta = span.first; // updating
                lastR2AlignEnd = span.second;                 
            }
            else ;            
            // update current as last
            lastAlignSta = span.first;
            lastAlignEnd = span.second;
        }
        // summary last cluster: first_alignment_start all_alignment_cov_len	base_cov	read_num
        assert(lastAlignEnd>firstAlignSta);
        moleculeLen = lastAlignEnd-firstAlignSta+1;
        if(moleculeLen==0) moleculeLen = 1; // caution: avoid potential dividing 0
        moutfp << lastAlignSta     << "\t" 
               << moleculeLen      << "\t" 
               << setprecision(2)  << fixed << (float)covleng/moleculeLen << "\t" 
               << readNum          << "\tU\t" 
               << lastAlignEnd     << "\t"
               << alignedPos.str() << "\t"
               << rpairids.str()   << "\t" 
               << readids.str()    << "\t"
               << rcigars.str()    << endl;
        if(moleculeLen >= minMoleSize)
        {
            // molecule cov
            std::stringstream mckey;
            mckey.str("");
            mckey << setprecision(2) << fixed << (float)covleng/moleculeLen;
            map<string, unsigned long>::iterator mcitr = (*molecule_cov).find(mckey.str());
            if(mcitr == (*molecule_cov).end())
            {
                (*molecule_cov).insert(std::pair<string, unsigned long>(mckey.str(), 1));
            }
            else
            {
                (*mcitr).second += 1;
            }
            // molecule len
            int mlkey = (int)round((double)moleculeLen/1000); 
            map<int, unsigned long>::iterator mlitr = (*molecule_len).find(mlkey);
            if(mlitr == (*molecule_len).end())
            {
                (*molecule_len).insert(std::pair<int, unsigned long>(mlkey, 1));
                (*molecule_len_sum).insert(std::pair<int, unsigned long>(mlkey, moleculeLen));            
            }
            else
            {
                (*mlitr).second += 1;
                (*molecule_len_sum)[mlkey] += moleculeLen;            
            }
            // molecules per barcode
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
            // molecule read number
            map<int,    unsigned long>::iterator readitr;
            readitr = (*molecule_read_num).find(readNum);
            if(readitr == (*molecule_read_num).end())
            {
                (*molecule_read_num).insert(std::pair<int, unsigned long>(readNum, 1));
                (*molecule_read_num_sum).insert(std::pair<int, unsigned long>(readNum, lastAlignEnd-lastAlignSta+1));                
            }
            else
            {
                (*molecule_read_num)[readNum] += 1;
                (*molecule_read_num_sum)[readNum] += lastAlignEnd-lastAlignSta+1;
            }
            // molecule read distances:
            map<int, unsigned long>::iterator imrditr;
            map<int, unsigned long>::iterator imrditr_end;
            map<int, unsigned long>::iterator imrditr_final;
            imrditr     = im_read_interval.begin();
            imrditr_end = im_read_interval.end();
            while(imrditr != imrditr_end)
            {
                //if(abs((*imrditr).first) > 500) // note: with 0.1x of 50kb molecule, dist between 18 readpairs is 2777bp
                if(true)
                {
                    imrditr_final = (*molecule_read_interval).find((*imrditr).first);
                    if(imrditr_final == (*molecule_read_interval).end())
                    {
                        (*molecule_read_interval).insert(std::pair<int, unsigned long>((*imrditr).first, (*imrditr).second));
                    }
                    else
                    {
                        (*imrditr_final).second += (*imrditr).second;
                    }
                }
                imrditr ++;
            }
        }
        //
        im_read_interval.clear();    
        readNum = 0;
        //
        alignedPos.str("");
        // read id and order of pair info
        readids.str("");
        rcigars.str("");
        rpairids.str("");        
        mitr = clitr;     
    }
    cout << "   Info: " << lastChr << " has " << alignNum << " alignments analyzed." << endl;
    return true;
}

bool write_molecule_stat(map<string, unsigned long> molecule_cov, 
                         map<int,    unsigned long> molecule_len, 
                         map<int,    unsigned long> molecule_len_sum, 
                         map<int,    unsigned long> molecule_read_num,
                         map<int,    unsigned long> molecule_read_num_sum,
                         map<int,    unsigned long> molecule_read_interval,
                         map<string, int>           barcode_molecule_cnt,
                         string out_prefix,
                         int minMoleSize)
{
    // 1. molecule base coverage
    std::stringstream minMsize;
    minMsize.str("");
    minMsize << minMoleSize;
    string outMCovFile = out_prefix + "_min" + minMsize.str() + "bp_moleCov_stat.txt";
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
    string outMLenFile = out_prefix + "_min" + minMsize.str() + "bp_moleLen_stat.txt";    
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
    string outMpBHistFile = out_prefix + "_min" + minMsize.str() + "bp_moleNumPerBarc_stat.txt";        
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
    string outReadNumFile = out_prefix + "_min" + minMsize.str() + "bp_readNum_stat.txt";    
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
    // 5. distance between reads along chr#barcode.
    string outReadDistFile = out_prefix + "_min" + minMsize.str() + "bp_readDist_stat.txt";    
    ofstream ofp5;
    ofp5.open(outReadDistFile.c_str(), ios::out);
    if(!ofp5.good())
    {
        cout << "   Error: cannot open " << outReadDistFile << " to write data. " << endl;
        return false;
    }     
    map<int,    unsigned long>::iterator disitr;
    map<int,    unsigned long>::iterator disitr_end; 
    disitr     = molecule_read_interval.begin();
    disitr_end = molecule_read_interval.end();
    ofp5 << "#read_distance\tnumber_of_cases_with_that_read_number" << endl;
    while(disitr != disitr_end)
    {
        ofp5 << (*disitr).first << "\t" << (*disitr).second << endl;
        disitr ++;
    }
    ofp5.close();
    return true;
}

int decipher_cigar(string cigar)
{
    // CIGAR alphabet: M I D N S H P; * = X: 23S23M1D8M1I16M1D29M
    char *cstr = (char*)cigar.c_str();
    string numstr("");
    int covlen = 0;
    for (int i=0; i<cigar.size(); i++)
    {
        if(cstr[i]>='0' && cstr[i]<='9')
        {
            numstr += cstr[i];
        }
        else
        if(cstr[i] == 'M' || cstr[i] == 'I'|| cstr[i] == 'D'|| cstr[i] == 'N'|| 
           cstr[i] == 'S' || cstr[i] == 'H'|| cstr[i] == 'P')
        {
            if(cstr[i] == 'M' || cstr[i] == 'D')
            {
                covlen += atoi(numstr.c_str());
            }
            numstr.clear();
        }
    }    
    return covlen;
}


