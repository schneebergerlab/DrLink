/* this function, given a known marker file (labeled with parental info) and a vcf file from longranger wgs
   output the distribution of breadpoints of recombination.
   
   this is different from previous version on molecule analysis with marker; it starts with molecules in recovered from bam.
   2018-01-01
   
   this checks alignments on selecting markers in recombis step.
*/
#include                   <stdio.h>
#include                  <string.h>
#include                    <string>
#include                  <stdlib.h>
#include                  <iostream>
#include                    <time.h>
#include                  <assert.h>
#include                   <iomanip>
#include            "split_string.h"
#include             "convert_bam.h"
#include        "analyze_molecule.h"
#include       "detect_recombsite.h"
#include       "find_intersection.h"
#include         "sample_molecule.h"
#include      "sample_molecule_v2.h"
#include    "separate_kb_molecule.h"

using namespace std;

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cout << endl;
        cout << "   Given linked-read sequencing of recombinants, this tool identifies meiotic recombinations." << endl;
        const char *buildString = __DATE__ ", " __TIME__ "";
        cout << "   (version 1.0 - compiled on " << buildString << ")"                   << endl  << endl;
        cout << "   Usage: DrLink subprogram [options] "                                          << endl;
        cout << "\n   Subprograms: "                                                              << endl;
        cout << "   \tpreprocess \t(Usage: please check options with DrLink preprocess)"          << endl;  
        cout << "   \tmolecule   \t(Usage: please check options with DrLink molecule)"            << endl;
        cout << "   \trecombis   \t(Usage: please check options with DrLink recombis)"            << endl;
        cout << "   \tintersec   \t(Usage: please check options with DrLink intersec)"    << endl << endl;        
        exit(1);
    }
    double startT= clock();
    string subprogram = (string)argv[1];
    // step 1. select info from bam with samtools view
    if(subprogram.compare("preprocess")==0 && !convertBam(argc, argv))
    {
        if(argc < 4)
        {
            cout << "\n   Info: not enough inputs. "          << endl << endl;
        }
        else
        {
            cout << "   Error: no bam info coverted - please check you bam file with BX? " << endl << endl;
        }
        return false;
    }
    // step 2. analyze molecules with the selected information
    if(subprogram.compare("molecule")==0 && !analyze_molecule(argc, argv))
    {
        if(argc < 6) 
        {
            cout << "\n   Info: not enough inputs. "          << endl << endl;
        }
        else
        {
            cout << "   Error: failed in molecule analysis. " << endl << endl;
        }
        return false;
    }   
    // step 3. detect recombination sites with vcf and known markers
    if(subprogram.compare("recombis")==0 && !detect_recombsite(argc, argv))
    {
        if(argc < 13) 
        {
            cout << "\n   Info: not enough inputs. "          << endl << endl;
        }
        else
        {
            cout << "   Error: failed in recombination site detection. " << endl << endl;
        }
        return false;
    }     
    // func 4. find intersection of two sets of recombination site predictions
    if(subprogram.compare("intersec")==0 && !find_intersection(argc, argv))
    {
        if(argc < 7) 
        {
            cout << "\n   Info: not enough inputs. "          << endl << endl;
        }
        else
        {
            cout << "   Error: failed in finding intersection of the two predictions. " << endl << endl;
        }
        return false;
    }
    // func 5. sample molecules according to molecule size distribution (for >=2 pools)
    if(subprogram.compare("resample")==0 && !sample_molecule(argc, argv))
    {
        if(argc < 7) 
        {
            cout << "\n   Info: not enough inputs. "          << endl << endl;
        }
        else
        {
            cout << "   Error: failed in sampling molecules. " << endl << endl;
        }
        return false;
    }    
    // func 6. sample molecules - with both size and base coverage
    if(subprogram.compare("separate_kb")==0 && !separate_kb_molecule(argc, argv))
    {
        if(argc < 5) 
        {
            cout << "\n   Info: not enough inputs. "          << endl << endl;
        }
        else
        {
            cout << "   Error: failed in sampling molecules. " << endl << endl;
        }
        return false;
    }  
    // func 7. sample molecules according to base coverage distribution for 1kb,2kb,...
    if(subprogram.compare("resample_v2")==0 && !sample_molecule_v2(argc, argv))
    {
        if(argc < 8)
        {
            cout << "\n   Info: not enough inputs. "          << endl << endl;
        }
        else
        {
            cout << "   Error: failed in sampling molecules. " << endl << endl;
        }
        return false;
    }        
    //
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    return 0;
}
