/* this function select info from a bam file, which is piped from samtools 

   new on 2018-06-07:
   1. remove discortly aligned read pairs (inter-chromosome ones).
   
*/

#include  <iostream>
#include    <vector>
#include    <string>
#include  <stdlib.h>
#include    <time.h>  /* clock_t, clock, CLOCKS_PER_SEC */
#include  "./gzlib/gzstream.h"

#include "split_string.h"

bool convertBam(int argc, char* argv[])
{
    if(argc < 4)
    {
        cout << "\n   DrLink preprocess\t"
             << "Usage: samtools view longranger_wgs.bam | DrLink preprocess - out_prefix_str " << endl;  
        return false;
    }
    clock_t tbeg;
    tbeg = clock();
    bool verbose = true; // to add future option    
    //
    cout << "   Info: converting streaming bam info started..." << endl;
    // prepare output gzipped file
    string ofilename = (string)argv[3] + "_bamInfo_trashme.gz";
    ogzstream  ofp;
    ofp.open(ofilename.c_str());
    if ( ! ofp.good()) {
        std::cerr << "ERROR: Opening file `" << ofilename << "' failed.\n";
	return 1;
    }
    ofp << "#chr\tpos\tbamCIGAR\tbarcode\tlrMI\tread-order:bamFLAG\treversed-align\tread-id\tpnext" << endl;
    // pnext: where its pair is aligned: must be the same-chr-position - only concordant alignments kept.
    //
    unsigned long numraw  = 0;
    unsigned long numcov  = 0;
    unsigned long numNoBX = 0;
    unsigned long numNoCG = 0; // no cigar string    
    unsigned long numNoMI = 0;
    unsigned long numNoMA = 0; // number of reads showing multiple alignment: only one kept.    
    unsigned long numDiscordant = 0; // read pairs aligned not in the same chr-region.
    string exampleNoBX("");
    string exampleNoCG("");    
    string exampleNoMI("");
    string exampleNoMA("");
    string exampleDiscordant("");
    bool out = false;
    std::string line;
    while (std::getline(std::cin, line)) 
    {
        if(line.compare("quit")==0 || line.compare("q")==0 || line.compare("exit")==0) break;
        if(line.size()==0 || line[0]=='#') continue;
        numraw ++;
        if(numraw%10000000 == 0)
        {
            cout << "   info: " << numraw << "th aligm..." << endl;
        }
        //
        vector<string> lineinfo = split_string(line, '\t');
        
        // cigar string
        if(lineinfo[5].compare("*") == 0) 
        {
            if(exampleNoCG.size()==0)
            {
                exampleNoCG = line;
                cout << "   Warning: there are alignments without explicit CIGAR string, skipped., e.g.: " 
                     << line 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;
            }
            numNoCG ++;
            continue;
        }
        // discordant read pairs
        if(lineinfo[6].compare("=") != 0)
        {
            // case 1: not on the same chr
            if(exampleDiscordant.size()==0)
            {
                exampleDiscordant = line;
                cout << "   Warning: there are discordant alignments, skipped., e.g.: " 
                     << line 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;
            }  
            numDiscordant ++;
            continue;                      
        }
        else
        {
            // case 2: on the same read, but unexpected template length 
            string tmplen;
            if(lineinfo[8].find("-") != std::string::npos)
            {
                tmplen = lineinfo[8].substr(1);
            }
            else
            {
                tmplen = lineinfo[8].substr(0);
            }
            unsigned long templateLen = strtoul(tmplen.c_str(), NULL, 0);
            if(templateLen==0 || templateLen>2000) // caution: fragment size
            {
                if(exampleDiscordant.size()==0)
                {
                    exampleDiscordant = line;
                    cout << "   Warning: there are discordant alignments, skipped., e.g.: " 
                         << line 
                         << endl;
                    cout << "   Info: you will get a total number of such alignments in the end of program. " 
                         << endl;
                }
                numDiscordant ++;
                continue;    
            }
        }
        //
        /* BX: Chromium barcode sequence that is error-corrected and confirmed against 
           a list of known-good barcode sequences. Use this for analysis, e.g., '\tBX:Z:TTTGTGTGTATGACTC-1' */        
        size_t pos1  = line.find("\tBX:");
        if(pos1 == std::string::npos) 
        {
            if(exampleNoBX.size()==0)
            {
                exampleNoBX = line;
                cout << "   Warning: there are alignments without BX info, skipped., e.g.: " 
                     << line 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;                
            }
            numNoBX ++;
            continue;
        }
        pos1 += 1; // from 'BX:Z:...'
        size_t pos2  = line.find("\t", pos1);
        if(pos2 == std::string::npos)
        {
            pos2 = line.size();
        }        
        string rawbx = line.substr(pos1, pos2-pos1);
        //
        /* MI: Global molecule identifier for molecule that generated this read,
           e.g.,  '\tMI:i:29' */
        size_t pos3  = line.find("\tMI:");
        bool foundmi = true;
        if(pos3 == std::string::npos)
        {
            if(exampleNoMI.size()==0)
            {
                exampleNoMI = line;
                cout << "   Warning: there are alignments without MI info, but kept, e.g.: " 
                     << line 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;                
            }
            numNoMI ++;
            foundmi = false;
        }
        pos3 += 1; // from 'MI:i:...'
        size_t pos4;
        string rawmi("-1");
        if(foundmi)
        {
            pos4 = line.find("\t", pos3);
            if(pos4 == std::string::npos)
            {
                pos4 = line.size();
            }
            rawmi = line.substr(pos3+5, pos4-pos3-5); // without "MI:i:"
        }
        //
        // check if it is read 1 or read 2
        int hexflag = strtol(lineinfo[1].c_str(), NULL, 0);
        if(hexflag > 255) 
        {
            numNoMA ++;
            if(exampleNoMA.size()==0)
            {
                exampleNoMA = line;
                cout << "   Warning: there are alignments being secondary/supplementary or not passing filters, skipped, e.g.: " 
                     << line 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;                
            }
            continue;
        }

        string readflag = "U";
        
        if((hexflag & 0x40) == 0x40 && (hexflag & 0x80) == 0)
        {
            readflag = "R1";
        }
        else
        if((hexflag & 0x40) == 0    && (hexflag & 0x80) == 0x80)
        {
            readflag = "R2";
        }
        else
        {
            readflag = "U0";
            cout << "   hexflag=" << hexflag << endl;
            cout << "   hexflag & 0x40 == 0x40= " << (hexflag & 0x40) << endl;
            cout << "   hexflag & 0x80 == 0   = " << (hexflag & 0x80) << endl;            
        }
        string reversed = "nrc";
        if((hexflag & 0x10) == 0x10) reversed = "rc";
        // chr	pos	CIGAR	barcode	MI	read-ordering:R1/R2	reversed	read-id	pnext
        // caution!!! code below also works for reversed-cases of alignment!
        // in bam/sam, it is always pos + CIGAR"M" ==> spanning of a read!!!
        ofp << lineinfo[2]                     << "\t" 
            << lineinfo[3]                     << "\t" 
            << lineinfo[5]                     << "\t" 
            << rawbx.substr(5, rawbx.size()-7) << "\t"
            << rawmi                           << "\t"
            << readflag  << ":" << lineinfo[1] << "\t" 
            << reversed                        << "\t"
            << lineinfo[0]                     << "\t"
            << lineinfo[7]                     << endl; // with read id and pnext
        //
        numcov ++;            
        if(!out) out = true;
        line.clear();
    }    
    //
    ofp.close();

    if(numNoBX > 0)
    {
        cout << "   Warning: there are "<< numNoBX 
             << " alignments without BX info, skipped. "              
             << endl;
    }
    if(numNoCG > 0)
    {
        cout << "   Warning: there are "<< numNoCG 
             << " alignments without explicit CIAGR info, skipped."   
             << endl;     
    }    
    if(numNoMI > 0)
    {
        cout << "   Warning: there are "<< numNoMI 
             << " alignments without MI info, but kept."              
             << endl;     
    }
    if(numNoMA > 0)
    {
        cout << "   Warning: there are "<< numNoMA 
             << " alignments being secondary/supplementary alignment"
             << " or not passing filters, skipped." 
             << endl;     
    }  
    if(numDiscordant > 0)
    {
        cout << "   Warning: there are "<< numDiscordant 
             << " alignments being discordant (pairs not on same chr or distance>2kb), skipped." 
             << endl;         
    }  
    cout << "   Info: in total " 
         << numraw  
         << " aligment lines, of which " 
         << numcov               
         << " with barcode info coverted." 
         <<  endl;
    cout << "   Info: time on converting streaming bam info: " 
         << (float)(clock()-tbeg)/CLOCKS_PER_SEC 
         << " seconds.\n" 
         << endl;
    if(!out) return false;
    else return true;
}
