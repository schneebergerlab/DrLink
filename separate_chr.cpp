#include        "./gzlib/gzstream.h"
#include                    <zlib.h>
#include                   <stdio.h>
#include                  <string.h>
#include                    <string>
#include                       <map>
#include                    <vector>
#include                  <stdlib.h>
#include                  <iostream>
//#include            "split_string.h"
using namespace std;
bool gen_random(string* str, const int len);
//
bool separate_chr(string filename, map<string, string>* visitedChr, string* tmpflag)
{
    /* this function separates 
    
       xxx_molecule_table_trashme.txt.gz from DrLink molecule or 
       phased_variants.vcf.gz            from longranger
       
       into chr-sepecific sub-files for later processing.
       
       Input:
       
       filename  - name of input gzipped file.
       moleORvcf - if "moleTable", the file is treated as molecule table file; if "vcf", vcf file.
       
       return variable:
       
       visitedChr- map<chr, chrfilename>
    */
    // files need to be gzipped
    if(filename.size() < 4)
    {
        cout << "   Error: unexpted file name " << filename << ". Must be in form of xxx.gz ." << endl;
        return false;
    }
    if(filename.substr(filename.size()-3).compare(".gz") != 0)
    {
        cout << "   Error: only file xxx.gz for molecule table or vcf accepted. " << endl;
        return false;
    }
    //
    string tmpfileflag= *tmpflag;
    // handle of utput file
    ogzstream outfp;
    string lastChr("");
    // read and write
    igzstream  fp;
    fp.open(filename.c_str());
    if(!fp.good())
    {
        cout << "   Error: cannot open file " << filename << "; exited. " << endl;
        return false;
    }
    while(fp.good())
    {
        string line("");
        getline(fp, line);
        if(line.size()==0 || line[0]=='#') continue;
        
        //vector<string> lineinfo = split_string(line, '\t');
        
        string chrinfo(""); // == lineinfo[0]
        size_t tabfirst = line.find("\t");
        chrinfo = line.substr(0, tabfirst);
        
        // find chr
        size_t ipos = chrinfo.find("#");     // if found - molecule table; if not - vcf file
        bool   vcf  = false;
        if(ipos == std::string::npos)
        {
            ipos = chrinfo.size();
            vcf  = true;
        }
        string chr;
        if(line.substr(0, 3).compare("Chr")==0 || 
           line.substr(0, 3).compare("chr")==0 || 
           line.substr(0, 3).compare("CHR")==0)
        {
            chr = chrinfo.substr(3, ipos-3); 
        }
        else
        {
            chr = chrinfo.substr(0, ipos); 
        }        
        // if same as last Chr, directly output
        if(lastChr.size()>0 && chr.compare(lastChr)==0)
        {
            if(line.substr(0, 3).compare("Chr")==0 || 
               line.substr(0, 3).compare("chr")==0 || 
               line.substr(0, 3).compare("CHR")==0)
            {
                outfp << line.substr(3) << endl;
            }
            else
            {
                outfp << line << endl;
            }
            continue;
        }
        // otherwise, close file for last Chr; open new file for current chr if not exist
        if((*visitedChr).find(chr) == (*visitedChr).end())
        {
            string chrfilename("");
            if(vcf == true)
            {
                int t=0;
                chrfilename = tmpfileflag + "_subfile_chr_" + chr + "_longranger_variants.vcf.gz";
                igzstream iitfp;
                iitfp.open(chrfilename.c_str(), ios::in);
                while(t<10000 && iitfp.good())
                {
                    iitfp.close();
                    tmpfileflag.clear();
                    gen_random(&tmpfileflag, 6);
                    chrfilename = tmpfileflag + "_subfile_chr_" + chr + "_longranger_variants.vcf.gz";
                    iitfp.open(chrfilename.c_str(), ios::in);
                    if(iitfp.good())
                    cout << "   Warning: tmp file " << chrfilename << " exists! Trying another!" << endl;
                }
            }
            else
            {
                int t=0;
                chrfilename = tmpfileflag + "_subfile_chr_" + chr + "_DrLink_molecules.txt.gz";
                igzstream iitfp;
                iitfp.open(chrfilename.c_str(), ios::in);
                while(t<10000 && iitfp.good())
                {
                    iitfp.close();
                    tmpfileflag.clear();
                    gen_random(&tmpfileflag, 6);
                    chrfilename = tmpfileflag + "_subfile_chr_" + chr + "_DrLink_molecules.txt.gz";
                    iitfp.open(chrfilename.c_str(), ios::in);
                    if(iitfp.good())
                    cout << "   Warning: tmp file " << chrfilename << " exists! Trying another!" << endl;
                }                
            }
            if(outfp.good())
            {
                outfp.close();
            }
            outfp.open(chrfilename.c_str(), ios::out);
            cout << "   Info: new subfile " << chrfilename << " created." << endl;
            // output
            if(line.substr(0, 3).compare("Chr")==0 || 
               line.substr(0, 3).compare("chr")==0 || 
               line.substr(0, 3).compare("CHR")==0)
            {
                outfp << line.substr(3) << endl;
            }
            else
            {
                outfp << line << endl;
            }
            // collect chr and filename
            (*visitedChr).insert(std::pair<string, string>(chr, chrfilename));
        }
        else
        {
            string chrfilename = (*visitedChr)[chr];
            if(outfp.good())
            {
                outfp.close();
            }
            outfp.open(chrfilename.c_str(), ios::out | ios::app);
            cout << "   Info: subfile " << chrfilename << " re-opened." << endl;                 
            // output 
            if(line.substr(0, 3).compare("Chr")==0 || 
               line.substr(0, 3).compare("chr")==0 || 
               line.substr(0, 3).compare("CHR")==0)
            {
                outfp << line.substr(3) << endl;
            }
            else
            {
                outfp << line << endl;
            }
        }
        // update lastChr as current chr for next line
        lastChr = chr;
    }
    fp.close();
    outfp.close();
    
    (*tmpflag) = tmpfileflag;
    return true;
}

bool gen_random(string* str, const int len) {
    // from stackoverflow
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    for (int i = 0; i < len; ++i) {
        (*str) += alphanum[rand() % (sizeof(alphanum) - 1)];
    }
    return true;
}
