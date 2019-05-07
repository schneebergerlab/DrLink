/* This function aims to remove false positives resulting from repetitiveness in genomes.

   Given two sets of CO predictions in the same read set based on parental genomes, 
   find the intersection of them.
   
   Four files will be output:
   1. commonly found: pa coordinate
   2. commonly found: pa + pb coordinate
   3. specific       CO in   pa
   4. specific       CO in   pb 

   Hequan Sun
   Date: 2018-08-28 15:34
*/
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
#include                  "globals.h"
#include             "split_string.h"
#include       "check_allele_order.h"
#include            "insert_raw_bp.h"
#include "postprocess_bpPrediction.h"
#include                "output_co.h"

using namespace std;
//
bool get_CO(string cofile, multimap<string, BPINFO>* aco);
bool get_CO_reads(string coreadfile, map<string, int>* reads);
//
bool find_intersection(int argc, char* argv[])
{
    std::stringstream usage;
    usage.str("");
    if(argc < 7)
    {
        usage << endl;
        usage << "   Given two sets of CO predictions using parents, this function finds out the intersection." << endl;
        const char *buildString = __DATE__", " __TIME__;
        usage << "   (compiled on " << buildString << ")"                                               << endl << endl;
        usage << "   Usage: DrLink intersec pa-CO.txt pb-CO.txt pa-co-reads pb-co-reads outfolder "     << endl << endl;
        usage << "\twhere "                                                         << endl;
        usage << "\tpa-CO.txt:    CO prediction using parent a/Col as reference"    << endl;
        usage << "\tpb-CO.txt:    CO prediction using parent b/Ler as reference"    << endl;
        usage << "\tpa-co-reads:  CO-read path given  by recombis for parent a/Col" << endl;
        usage << "\tpb-co-reads:  CO-read path given  by recombis for parent b/Ler" << endl;        
        usage << "\toutfolder:    output folder"                                    << endl;
        cout << usage.str() << endl;
        return false;
    }
    double startT= clock();   
    string outfolder = (string)argv[6]+"_intersection";
    if(outfolder.size()==0)   
    outfolder = "enjo_intersection";
    // create an intermediate folder for collecting details about a CO-molecule
    DIR* dir = opendir(outfolder.c_str());
    if (dir)
    {
        /* Directory exists. */
        closedir(dir);
    }
    else if (ENOENT == errno)
    {
        /* Directory does not exist. */
        const int dir_err = mkdir(outfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (dir_err == -1)
        {
            cout << "   Error: cannot create directory " << outfolder << endl;
            return false;
        }        
    }
    else ;
    // 1.read parent-referenced COs   
    string acofile     = argv[2];
    string bcofile     = argv[3];
    string acoreadpath = argv[4];
    string bcoreadpath = argv[5];
    multimap<string, BPINFO> aco; // <chr#barcode, BPINFO>
    if(!get_CO(acofile, &aco))
    {
       return false;
    }
    multimap<string, BPINFO> bco;
    if(!get_CO(bcofile, &bco))
    {
       return false;
    }    
    // 2. find common COs in both findings
    multimap<string, BPINFO>::iterator cbitr;
    multimap<string, BPINFO>::iterator cbitr_end;
    cbitr     = aco.begin();
    cbitr_end = aco.end();
    while(cbitr != cbitr_end)
    {
        string key   = (*cbitr).first;
        BPINFO tmpbp = (*cbitr).second;
        // get "a" reads: filename like: path/to/1_1132082_1145010_12_CGGCTAGTCGCCGTGA.txt
        vector<string> keyinfo = split_string(key, '#');
        std::stringstream acoreadfile;
        acoreadfile.str("");
        acoreadfile << acoreadpath  << "/" 
                    << keyinfo[0]   << "_" 
                    << tmpbp.leftp  << "_"
                    << tmpbp.rightp << "_"
                    << tmpbp.order  << "_"
                    << keyinfo[1]   << ".txt\0";
        map<string, int> aReads;
        if(!get_CO_reads(acoreadfile.str(), &aReads))
        {
            return false;
        }
        //
        std::pair <std::multimap<string, BPINFO>::iterator, std::multimap<string, BPINFO>::iterator> ret;
        ret = bco.equal_range(key);
        multimap<string, BPINFO>::iterator tmpitr; 
        for(tmpitr = ret.first; tmpitr != ret.second; tmpitr ++)
        {
            // get "b" reads
            string         key2     = (*tmpitr).first;
            BPINFO         tmpbp2   = (*tmpitr).second;
            vector<string> keyinfo2 = split_string(key2, '#');
            std::stringstream bcoreadfile;
            bcoreadfile.str("");
            bcoreadfile << bcoreadpath   << "/" 
                        << keyinfo2[0]   << "_" 
                        << tmpbp2.leftp  << "_"
                        << tmpbp2.rightp << "_"
                        << tmpbp2.order  << "_"
                        << keyinfo2[1]   << ".txt\0";
            map<string, int> bReads;
            if(!get_CO_reads(bcoreadfile.str(), &bReads))
            {
                return false;
            }
            // check if "a" reads match "b" reads
            map<string, int>::iterator ritr;
            map<string, int>::iterator ritr_end;
            ritr      = aReads.begin();
            ritr_end  = aReads.end();
            int found = 0;
            while(ritr != ritr_end)
            {
                if(bReads.find((*ritr).first) != bReads.end())
                {
                    found ++;
                }
                ritr ++;
            }
            if(found>=3)
            {
                //
                if( (*cbitr).second.ovline.size() > 0)
                {
                    (*cbitr).second.ovline += "\n";
                }
                (*cbitr).second.ovline += tmpbp2.bpline;
                //
                if( (*tmpitr).second.ovline.size() > 0)
                {
                    (*tmpitr).second.ovline += "\n";
                }
                (*tmpitr).second.ovline += tmpbp.bpline;
            }
        }
        cbitr ++;
    }    
    // 3. output files
    // collect "a/Col" CO predictions for sorting and overlapping
    multimap<string, multimap<unsigned long, BPINFO> > commonbpraw_acoord;
    // file with commonly found co in both predictions => both coordinates for checking.
    string ofilename2 = outfolder + "/ref_CL_common_bp.txt";
    ofstream ofp2;
    ofp2.open(ofilename2.c_str(), ios::out);
    if(!ofp2.is_open())
    {
        cout << "   Error: cannot open file " << ofilename2 << endl;
        return false;
    }  
    // file with co specific to a/Col
    string ofilename3 = outfolder + "/ref_C__specific_bp.txt";
    ofstream ofp3;
    ofp3.open(ofilename3.c_str(), ios::out);
    if(!ofp3.is_open())
    {
        cout << "   Error: cannot open file " << ofilename3 << endl;
        return false;
    }  
    // collecting
    cbitr     = aco.begin();
    cbitr_end = aco.end();
    while(cbitr != cbitr_end)
    {
        string key   = (*cbitr).first;
        BPINFO tmpbp = (*cbitr).second;
        if(tmpbp.ovline.size() > 0)
        {
            insert_raw_bp(&commonbpraw_acoord, 
                          tmpbp.chr, 
                          tmpbp.leftp, 
                          tmpbp.rightp, 
                          tmpbp.score,
                          tmpbp.order, 
                          tmpbp.bpline);            
            ofp2 << tmpbp.bpline   << "\tref.C\n";   // "a/Col"
            ofp2 << tmpbp.ovline   << "\tref.L\n\n"; // "b/Ler"
        }   
        else
        {
            ofp3 << tmpbp.bpline   << endl;
        }     
        cbitr ++;
    }
    // file with co specific to b/Ler
    string ofilename4 = outfolder + "/ref_L__specific_bp.txt";
    ofstream ofp4;
    ofp4.open(ofilename4.c_str(), ios::out);
    if(!ofp4.is_open())
    {
        cout << "   Error: cannot open file " << ofilename4 << endl;
        return false;
    } 
    cbitr     = bco.begin();
    cbitr_end = bco.end();
    while(cbitr != cbitr_end)
    {
        string key   = (*cbitr).first;
        BPINFO tmpbp = (*cbitr).second;
        if(tmpbp.ovline.size() <= 0)
        {
            ofp4 << tmpbp.bpline << endl;           
        } 
        cbitr ++;
    }        
    // close files
    ofp2.close();
    ofp3.close();
    ofp4.close();
    // raw non-overlapped bp - left-coordinate sorted raw break points
    string ofilename1 = outfolder + "/ref_C__common_bp_sorted.txt";
    ofstream ofp1;
    ofp1.open(ofilename1.c_str(), ios::out);
    if(!ofp1.is_open())
    {
        cout << "   Error: cannot open file " << ofilename1 << endl;
        return false;
    }    
    if(!output_co(commonbpraw_acoord, ofilename1))
    {
        return false;
    }
    // post-process predicated break points: final bp file
    string finalbpfile = outfolder +"/ref_C__common_bp_overlapped_final.txt";
    if(!postprocess_bpPrediction(commonbpraw_acoord, finalbpfile))
    {
        cout << "   Error: processing raw break points failed. " << endl;
        return false;
    }
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return true;
}
//
bool get_CO_reads(string coreadfile, map<string, int>* reads)
{
    ifstream ifp;
    ifp.open(coreadfile.c_str(), ios::in);
    if(!ifp.is_open())
    {
        cout << "   Error: cannot find file " << coreadfile << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        //00000000	0	1109895	1110045	2	J00137:120:HNWTCBBXX:2:1112:28970:19021	151M
        vector<string> lineinfo = split_string(line, '\t');
        (*reads).insert(std::pair<string, int>(lineinfo[5], 1));
    }
    ifp.close();
    return true;
}
//
bool get_CO(string cofile, multimap<string, BPINFO>* aco)
{
    cout << "   Info: reading raw bp info from file " << cofile << endl;
    ifstream ifp;
    ifp.open(cofile.c_str(), ios::in);
    if(!ifp.is_open())
    {
        cout << "   Error: cannot open file " << cofile << endl;
        return false;
    } 
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        vector<string> lineinfo = split_string(line, '\t');
        // at least with: 0.chr  1.left_pos  2.right_pos  3.score  4.left_alleles  5.right_alleles  6.barcode
        if(lineinfo.size() < 7) continue; 
        //
        BPINFO tmpbp;
        tmpbp.chr      = lineinfo[0];
        tmpbp.leftp    = strtoul(lineinfo[1].c_str(), NULL, 0);
        tmpbp.rightp   = strtoul(lineinfo[2].c_str(), NULL, 0);
        tmpbp.leftMax  = 0;
        tmpbp.rightMin = 0;
        tmpbp.bpline   = line;
        tmpbp.ovline   = "";
        tmpbp.score    = atof(lineinfo[3].c_str());
        tmpbp.order    = check_allele_order(lineinfo[4], lineinfo[5]);        
        // chr#bc#order
        string key = lineinfo[0] + "#" + lineinfo[6];              
        (*aco).insert(std::pair<string, BPINFO>(key, tmpbp));
    }
    ifp.close();  
    cout << "   Info: " << (*aco).size() << " bp info collected." << endl; 
    return true;
}
/*
                                        insert_raw_bp(&rawbpmap, 
                                                      chrbarcode[0], 
                                                      leftp, 
                                                      rightp, 
                                                      score,
                                                      aord, 
                                                      sstmp.str());

*/
