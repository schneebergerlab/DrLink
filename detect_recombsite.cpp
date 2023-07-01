/* this function detects recombination sites in meiosis, based on pooled linked read sequencing 

   input: 
      
      1. known markers from both parents
      2. variations called in pool     

   output:

      prediction on recombination sites.
      
   update: 2018-07-29 - check number of reads
*/
#include         "./gzlib/gzstream.h"
#include                     <zlib.h>
#include                    <stdio.h>
#include                   <string.h>
#include                     <string>
#include                        <map>
#include                     <vector>
#include                   <stdlib.h>
#include                    <fstream>
#include                   <iostream>
#include                    <sstream>
#include                     <time.h>
#include                   <assert.h>
#include                  <algorithm>
#include                    <iomanip>
#include                     <math.h>
#include                 <sys/stat.h>
#include                   <dirent.h>
#include                    <errno.h>
#include             "split_string.h"
#include             "separate_chr_v2.h"
#include            "insert_raw_bp.h"
#include       "check_allele_order.h"
#include "postprocess_bpPrediction.h"
#include                "output_co.h"
#include                  "globals.h"
//
struct iMARKER
{
    string         chr;
    unsigned long  pos;
    string        from;
    double          af; // observed alternative allele frequency: under high coverage, expect ~0.5, e.g., [0.3, 0.7]
    string         cov;
    bool       isGiven; // if the marker is in given list, set as true.
};
struct READINFO
{
    string readname;
    string mateid;
    unsigned long spansta;
    unsigned long spanend;
    string mallele; // marker1:allele1,\nmarker2:allele2
};
// --ivc 50 --iaf 0.2
unsigned long badBarcodeCov = 0;
int           minVarCov     = 50;    // minimum base coverage    of variant site (refcov+altcov) in vcf to be used as markers
int           maxVarCov     = 10000; // maximum base coverage    of variant site (refcov+altcov) in vcf to be used as markers
double        minVarAF      = 0.2;   // minimum allele frequency of variant site (refcov+altcov) in vcf to be used as markers
double        maxVarAF      = 1.0;   // maximum allele frequency of variant site (refcov+altcov) in vcf to be used as markers
int           minRead       = 2;     // minimum number of uniq read ids on each side of CO, excluding the one on CO
bool          rcm           = false;
multimap<string, unsigned long> HeteroBarcode;
string        tmpflagstr("tmptmp");
unsigned long total_check   = 0;
//
bool collect_marker_info(const char*           file, 
                         string                parent_id, 
                         map<string, string>*  mkr);
bool get_var_ref_barcode(string                formatInfo, 
                         string                barcodeString, 
                         vector<string>*       vbc,
                         vector<string>*       rbc,
                         string                this_chr,
                         string                this_pos,
                         multimap<string,      unsigned long>* HeteroBarcode);
bool get_break_pos_of_artSeq(string            artSeq, 
                         map<unsigned long,    string> artSeqPos, 
                         int                   min_marker,
                         unsigned long*        leftp, 
                         unsigned long*        rightp, 
                         int*                  bpPos,
                         double*               score);
bool get_moleStat_distribution(const char*     file, 
                         double                maxPercent, 
                         int*                  maxStat,
                         string                lenORreadnum);
bool get_bad_molecule(const char*              file, 
                         int                   maxReadNum,
                         int                   maxMoleLen,
                         multimap<string,      std::pair<unsigned long, unsigned long> >* badMolecule); // to remove 2018-01-01
bool read_chr_var(string                       vcffilename, 
                         map<string, string>   knownMarker, 
                         multimap<string,      iMARKER>* moleculeSet);
bool check_candidate_cosite(multimap<string,   std::pair<unsigned long, unsigned long> >* badMolecule,
                         string                mkey,
                         unsigned long         mstart,
                         unsigned long         mend);
bool build_artSeq(map<unsigned long, string>   alleleinfo, 
                         unsigned long         regionBegin, 
                         unsigned long         regionEnd, 
                         string*               artSeq,
                         string*               artSeqCheckingStr);                         
string check_allele_order(string               leftAlleles, 
                         string                rightAlleles);
bool get_recombis_options(int                  argc, 
                         char*                 argv[],
                         string*               filevcf,
                         string*               fileMoleTable,
                         string*               fileReadNumDistr,
                         string*               fileMoleLenDistr,
                         string*               filemark1,
                         string*               filemark2,
                         int*                  min_marker,
                         int*                  min_span,
                         double*               min_score,
                         string*               outPrefix,
                         int*                  max_moleculeLen,
                         int*                  max_readnum,
                         int*                  max_CORDist,
                         int*                  min_CORDist,
                         int*                  max_COSize,
                         map<string, string>* visitedVarChr,
                         map<string, string>* visitedMolChr);
bool check_co_read_density(unsigned long       co_start,
                         unsigned long         co_end,
                         string                read_align_info,
                         string                read_order_info,
                         string                read_id_info,
                         int*                  bp_per_read_in_co,
                         bool*                 bp_boundary_same_fragment,
                         int*                  min_RR_dist,
                         int*                  max_RR_dist,
                         int*                  n_readLeft,
                         int*                  n_readRight,
                         int*                  max_CORR_dist);
bool check_mole_local_base_cov(unsigned long   mole_start,
                         unsigned long         mole_end,
                         string                read_align_info,
                         unsigned long         win_size,
                         unsigned long         win_step,
                         double*               min_win_base_cov,
                         double*               max_win_base_cov,
                         double*               mean_win_base_cov,
                         double*               sd_win_base_cov);
bool output_real_molecule_len_distribution(
                         map<int, unsigned long> molecule_len_control, 
                         string outPrefix);
int get_chr_num_ii(string marker_file);
bool verbose = true; // TODO: get from option
// main
bool detect_recombsite(int argc, char* argv[])
{
    std::stringstream usage;
    usage.str("");
    if(argc < 14)
    {
        usage << endl;
        usage << "   Given necessary info, this function predicts meiotic recombination breakpoints."       << endl;
        const char *buildString = __DATE__", " __TIME__;
        usage << "   (special version compiled on " << buildString << ")"                           << endl << endl;
        usage << "   Usage: DrLink recombis options [default] "                                     << endl << endl;
        usage << "\t# Mandatory options: "                                                                  << endl;
        usage << "\t--var    STRING    snp variants                  in vcf.gz from longranger wgs  [NULL]" << endl;
        usage << "\t--mol    STRING    molecule meta info            in txt.gz from DrLink molecule [NULL]" << endl;
        usage << "\t--str    STRING    stat of molecule read numbers in txt    from DrLink molecule [NULL]" << endl;
        usage << "\t--stl    STRING    stat of molecule lengths      in txt    from DrLink molecule [NULL]" << endl;
        usage << "\t--spa    STRING    snp markers                   in txt    from paternal line   [NULL]" << endl;
        usage << "\t--sma    STRING    snp markers                   in txt    from maternal line   [NULL]" << endl;
        usage << "\t# Optional: "                                                                           << endl;
        usage << "\t--ivc    INT       min number of read coverage   for collecting a marker        [ 250]" << endl;
        usage << "\t--avc    INT       max number of read coverage   for collecting a marker        [ 400]" << endl;
        usage << "\t--iaf    DOUBLE    min alternative AF            for collecting a marker        [0.35]" << endl;
        usage << "\t--aaf    DOUBLE    max alternative AF            for collecting a marker        [0.65]" << endl;   
        usage << "\t--nur    INT       max molecule read number      for selecting good molecules   [NULL]" << endl;
        usage << "\t--nul    INT       max molecule length           for selecting good molecules   [NULL]" << endl;        
        usage << "\t--imn    INT       min number of allelic markers for reporting  a CO            [   3]" << endl;
        usage << "\t--imr    INT(bp)   min range  of allelic markers for reporting  a CO            [ 250]" << endl;        
        usage << "\t--ims    DOUBLE    min score                     for reporting  a CO            [ 0.8]" << endl;
        usage << "\t--icd    INT       min inter-read distance       for reporting  a CO            [ 500]" << endl;
        usage << "\t--acd    INT       max inter-read distance       for reporting  a CO            [5000]" << endl;
        usage << "\t--ird    INT       min number of side-reads      for reporting  a CO            [   2]" << endl;
      //usage << "\t--aco    INT       max interval size             for reporting  a CO            [2000]" << endl;
        usage << "\t--split  STRING    a flag refering to splitted --var and --mol                  [ off]" << endl;
        usage << "\t--rcm    bool      check if a marker is covered by a read                       [ off]" << endl;        
        usage << "\t--out    STRING    prefix                        for labeling output files      [enjo]" << endl;        
        usage                                                                                               << endl;
        usage << "\t*Note 1: one of the marker files can be null. "                                         << endl;
        usage << "\t*Note 2: if given --nur (--nul), --str (--strl) becomes invalid and ignored."           << endl;
        usage << "\t*Note 3: if given --split tmptmp and files found, then --var and --mol can be ignored. "<< endl;
        usage << "\t*Note 4: option --aco should be determined according to marker density\n\t(e.g., Athal: 300bp/marker), "
              << "and, read distribution shows a rough range of 0~5000 bp. " << endl;
        cout << usage.str() << endl;
        // --icd/acd is related to statistics of read coverage of 10x molecules
        // --ird: number of uniq read-ids beyond CO start and end, excluding the pair on CO.
        return false;
    }
    double startT= clock(); 
    // initialize variables
    string filevcf("");
    string fileMoleTable("");    
    string fileReadNumDistr("");
    string fileMoleLenDistr("");
    string filemark1("");
    string filemark2("");    
    int    min_marker = 3;    
    int    min_span   = 250;   // bp    
    double min_score  = 0.8;
    string outPrefix("enjo_");    
    int    maxReadNum = 0;     //22;     // if 0.1x per molecule: round(18*2*0.60)
    int    maxMoleLen = 0;     //54000;  // if 45kb per molecule: round(45000*2*0.60)
    int    maxCORDist = 5000;  //maximum inter-read distance in predicted CO intervals
    int    minCORDist = 500;
    int    maxCOSize  = 500000;//maximum size of interval to be predicted as a CO - no control now!!!
    map<string, string> visitedVarChr;
    map<string, string> visitedMolChr;
    // set variables from options
    if( !get_recombis_options(argc, 
                              argv,
                              &filevcf,
                              &fileMoleTable,
                              &fileReadNumDistr,
                              &fileMoleLenDistr,
                              &filemark1,
                              &filemark2,
                              &min_marker,
                              &min_span,
                              &min_score,
                              &outPrefix,
                              &maxMoleLen,
                              &maxReadNum,
                              &maxCORDist,
                              &minCORDist,
                              &maxCOSize,
                              &visitedVarChr,
                              &visitedMolChr) )
    {
        cout << "   Error: incorrect parameter settings. " << endl;
        return false;
    }
    multimap<string, std::pair<unsigned long, unsigned long> > badMolecule;
    // step m. read bad (longer or with larger read number) molecule info: decreasing fp, but increasing fn.
    // find max read number: set ratio of valley/peak*0.5 in "out_1000BC_90M_molecule_stat.pdf" after read peak
    if(maxReadNum == 0)
    if(!get_moleStat_distribution(fileReadNumDistr.c_str(), 0.20, &maxReadNum, "readNum"))
    {
        cout << "   Error: failed in reading file " << fileReadNumDistr << ". Exited." << endl;
        return false;
    }
    // find max molecule length: set ratio of valley/peak*0.5 in "out_1000BC_90M_molecule_stat.pdf" after molecule peak
    if(maxMoleLen == 0)
    if(!get_moleStat_distribution(fileMoleLenDistr.c_str(), 0.20, &maxMoleLen, "moleLen"))
    {
        cout << "   Error: failed in reading file " << fileMoleLenDistr << ". Exited." << endl;
        return false;
    }
    /* this is now done during reading of raw molecules
    if(!get_bad_molecule(fileMoleTable.c_str(), maxReadNum, maxMoleLen, &badMolecule))
    {
        cout << "   Error: failed in reading file " << fileMoleTable    << ". Exited." << endl;
        return false;            
    }
    */
    // step 0. initialized a variable for recording known marker info
    map<string, string> knownMarker;
    // step 1. read known marker info from parent 1
    if(!collect_marker_info(filemark1.c_str(), "1", &knownMarker))
    {
        cout << "   Error: cannot open marker file " << filemark1 << "; exited." << endl;
        exit(1);
    }
    // step 2. read known marker info from parent 2
    if(!collect_marker_info(filemark2.c_str(), "2", &knownMarker))
    {
        cout << "   Error: cannot open marker file " << filemark2 << "; exited." << endl;
        exit(1);
    }
    // step 3. separte molecule_table_trashme.txt.gz and phased_variants.vcf.gz into chromosome-specific files
    // seprate vcf-variant file: map<string, string> visitedVarChr
    if(visitedVarChr.size()==0 &&
       !separate_chr_v2(filevcf, &visitedVarChr, &tmpflagstr))
    {
        return false;
    }
    else
    {
        map<string, string>::iterator vchritr;
        map<string, string>::iterator vchritr_end;
        vchritr     = visitedVarChr.begin();
        vchritr_end = visitedVarChr.end();
        cout << "   Info: variants in vcf separated into chr-specific files: " << endl;
        while(vchritr != vchritr_end)
        {
            cout << "         " << (*vchritr).second << endl;
            vchritr ++;
        } 
    }
    // separate molecule file: map<string, string> visitedMolChr
    if(visitedMolChr.size()==0 &&
       !separate_chr_v2(fileMoleTable, &visitedMolChr, &tmpflagstr))
    {
        return false;
    }
    else
    {
        map<string, string>::iterator mchritr;
        map<string, string>::iterator mchritr_end;
        mchritr     = visitedMolChr.begin();
        mchritr_end = visitedMolChr.end();
        cout << "   Info: raw molecules separated into chr-specific files: " << endl;
        while(mchritr != mchritr_end)
        {
            cout << "         " << (*mchritr).second << endl;
            mchritr ++;
        }     
    }
    // step 5. prepare a file and a variable for collecting candidate break points
    std::stringstream ipre;
    ipre.str("");
    ipre << "_aML" << maxMoleLen
         << "_aRN" << maxReadNum
         << "_iSC" << min_score 
         << "_iMK" << min_marker
         << "_iSP" << min_span
         << "_iVC" << minVarCov
         << "_aVC" << maxVarCov
         << "_iAF" << minVarAF
         << "_aAF" << maxVarAF
         << "_ird" << minRead;
    // create an intermediate folder for collecting details about a CO-molecule
    DIR* dir = opendir((outPrefix+ipre.str()).c_str());
    if (dir)
    {
        /* Directory exists. */
        closedir(dir);
    }
    else if (ENOENT == errno)
    {
        /* Directory does not exist. */
        const int dir_err = mkdir((outPrefix+ipre.str()).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (dir_err == -1)
        {
            cout << "   Error: cannot create directory " << outPrefix+ipre.str() << endl;
            return false;
        }        
    }
    else ;
    // raw breakpoints file
    string rawbpfile = outPrefix + ipre.str() +"_raw_BP.txt";
    fstream rawfp;
    rawfp.open(rawbpfile.c_str(), ios::out);
    if(!rawfp.good())
    {
        cout << "   Error: cannot open file " << rawbpfile << " for collecting raw break point predictions. " << endl;
        return false;
    }
    rawfp << "#chr\t"
          << "left_pos\t"
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
          << "n_readLeft\t"
          << "n_readRight\t"
          << "max_CORR_dist\n";
    // files for collecting length of pure Col-/Ler-molecules: format: chr#barcode first-aligned-read mole_len artSeq
    string colMole = outPrefix + ipre.str() +"_col_mol_length.txt";
    string lerMole = outPrefix + ipre.str() +"_ler_mol_length.txt";
    fstream colMolefp;
    fstream lerMolefp;
    bool outLerCol = false;
    if(outLerCol)
    {
        colMolefp.open(colMole.c_str(), ios::out);
        lerMolefp.open(lerMole.c_str(), ios::out);
        if(!colMolefp.good())
        {
            cout << "   Error: cannot open file " << colMole << " for collecting molecule length info of Col. " << endl;
            return false;
        }
        if(!lerMolefp.good())
        {
            cout << "   Error: cannot open file " << lerMole << " for collecting molecule length info of Ler. " << endl;
            return false;
        }  
    }
    // file for collecting molecules passed all checkings on RN, ML and HET
    // for chi-square checking number of normal molecules in a window
    string goodmolefile = outPrefix + ipre.str() +"_pRNpMLpHET.txt.gz";
    ogzstream goodmolefp;
    goodmolefp.open(goodmolefile.c_str(), ios::out);
    if(!goodmolefp)
    {
        cout << "    Error: cannot open file " << goodmolefile << " for writing good molecule info. " << endl;
        return false;
    }
    goodmolefp << "#chr\tstart\tend\tRN" << endl;
    // a variable for collecting candidate break points
    multimap<string, multimap<unsigned long, BPINFO> > rawbpmap;
    //
    unsigned long numbad               = 0;
    unsigned long numbadML             = 0;
    unsigned long numbadRN             = 0;
    unsigned long numbadHET            = 0;
    unsigned long numbTotalMol         = 0;
    unsigned long numbTotalMolRead     = 0; // number of reads related to numbTotalMol
    unsigned long numbGoodTotalMol     = 0; // number of molecules passing checking on RN, ML and HET
    unsigned long numbGoodTotalMolRead = 0; // number of reads related to numbGoodTotalMol
    map<int, unsigned long> molecule_len_control; // length distribution of molecules passing all control
    //
    // step 6. analyze molecule for each chr, and allele info at markers according to vcf - 2018-01-01
    map<string, string>::iterator vchritr;
    map<string, string>::iterator vchritr_end;
    vchritr     = visitedVarChr.begin();
    vchritr_end = visitedVarChr.end();
    while(vchritr != vchritr_end)
    {
        string chr     = (*vchritr).first;
        cout << "   Info: finding breakpoints on chr " << chr << "..." << endl;
        // 
        cout << "      Extracting parental allele info at markers for molecules/barcodes... "    << endl;
        string vcffile = (*vchritr).second;
        multimap<string, iMARKER> moleculeSet; // <chr#bc, {marker-pos, from-parent} >
        if(!read_chr_var(vcffile, knownMarker, &moleculeSet))
        {
            cout << "      Error: read information at markers from " << vcffile << " failed. "   << endl;
            return false;
        }
        cout << "      Extracted parental allele info: " << moleculeSet.size() << "."            << endl;
        // find file of raw molecules for chr
        cout << "   Check: current chr " << chr << endl;
        assert(visitedMolChr.find(chr) != visitedMolChr.end());
        string molfile = visitedMolChr[chr];
        //
        igzstream molfp;
        molfp.open(molfile.c_str());
        if(!molfp.good())
        {
            cout << "   Error: cannot file of raw molecule for chr " <<  chr << ": " << molfile << endl;
            return false;
        }
        cout << "   Info: traversing the molecule table of chr " << chr << " to find breakpoints..." << endl;
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
            // filter away current molecule if with unexpected length or number of reads
            int    readnum = atoi(lineinfo[5].c_str());
            double molecov = atof(lineinfo[4].c_str());
            int    molelen = atoi(lineinfo[3].c_str()); 
            // pre-controlling on molecule length
            if(molelen<1000) continue;
            numbTotalMol ++;
            numbTotalMolRead += readnum;
            if(readnum>maxReadNum)
            {
                numbadRN ++;
                numbad ++;                
                continue;
            }
            else 
            if(molelen>maxMoleLen || 
               molelen<2*min_span)
            {
                // skip potential overlapping molecules, or short molecules
                numbadML ++;
                numbad ++;
                continue;
            }
            // find range of current molecule
            unsigned long msta = strtoul(lineinfo[1].c_str(), NULL, 0);
            unsigned long mend = strtoul(lineinfo[3].c_str(), NULL, 0) + msta -1;            
            // filter away current molecule if covering a heterozygous site
            std::pair <multimap<string, unsigned long>::iterator, 
                       multimap<string, unsigned long>::iterator> hetitr_region;
            hetitr_region = HeteroBarcode.equal_range(lineinfo[0]);
            multimap<string, unsigned long>::iterator hetitr;
            bool hetcov = false;
            for (hetitr=hetitr_region.first; hetitr!=hetitr_region.second; ++hetitr)
            {
                if((*hetitr).second>=msta && (*hetitr).second<=mend)
                {
                    cout << "   Info: molecule at " 
                         << lineinfo[0]      << ": " 
                         << lineinfo[1]      << "-" 
                         << lineinfo[2]      << " with het-bc at mkdr " 
                         << (*hetitr).second << " - skipped. " << endl;
                    numbad ++;
                    numbadHET ++;
                    hetcov = true;
                    break;
                }
            }
            if(hetcov == true) continue;
            // "real" moleclues passing checking on RN, ML and HET
            vector<string> chrbarcode = split_string(lineinfo[0], '#');          
            numbGoodTotalMol ++;
            numbGoodTotalMolRead += readnum;
            goodmolefp << chrbarcode[0] << "\t"
                       << lineinfo[1]   << "\t"
                       << lineinfo[2]   << "\t"
                       << lineinfo[5]   << endl;
            // collect length of a real molecule for calculating and checking expected CO molecule
            int mlkey = (int)round((double)molelen/1000); 
            map<int, unsigned long>::iterator mlitr = molecule_len_control.find(mlkey);
            if(mlitr == molecule_len_control.end())
            {
                molecule_len_control.insert(std::pair<int, unsigned long>(mlkey, 1));
            }
            else
            {
                (*mlitr).second += 1;
            }
            // new 2018-05-23
            //
            // for this molecule, get its alignment-covered positions for selecting markers
            //       8.all_reads_aligned_at  	15095089,15095230;15102519,15102668
            //       9.R1R2  		 	1,2
            //      10.read_id 		 	rid1,rid2
      	    /*    struct READINFO
		{
		    string readname;
		    string mateid;
		    unsigned long spansta;
		    unsigned long spanend;
		    string mallele; // marker1:allele1,\nmarker2:allele2
		};
	    */
            multimap<unsigned long, READINFO> molAlignMarkerAllele;
            
            map<unsigned long, string> molReadCoveredPos;
            vector<string> spaninfo  = split_string(lineinfo[8],  ';');            
            vector<string> mateinfo  = split_string(lineinfo[9],  ',');
            vector<string> nameinfo  = split_string(lineinfo[10], ',');
            vector<string> cigarinfo = split_string(lineinfo[11], ',');
            assert(spaninfo.size()==mateinfo.size() && mateinfo.size()==nameinfo.size());
            vector<string>::iterator rsp_itr     = spaninfo.begin();
            vector<string>::iterator rsp_itr_end = spaninfo.end();
            int geti = 0;
            bool readwarning = false;
            while(rsp_itr != rsp_itr_end)
            {
                string spanMidName = spaninfo[geti] + "\t" + mateinfo[geti] + "\t" + nameinfo[geti];
                //
                vector<string> indrsp = split_string(*rsp_itr, ',');
                unsigned long rstart  = strtoul(indrsp[0].c_str(), NULL, 0);
                unsigned long rend    = strtoul(indrsp[1].c_str(), NULL, 0);
                //
                if(rcm==true)
                for(unsigned long ii  = rstart; ii <= rend; ii ++)
                {
                    molReadCoveredPos.insert(std::pair<unsigned long, string>(ii, spanMidName));
                }
                //
                READINFO tmpallinfo;
                tmpallinfo.readname = nameinfo[geti]+"\t"+cigarinfo[geti];
                tmpallinfo.mateid   = mateinfo[geti];
                tmpallinfo.spansta  = rstart;
                tmpallinfo.spanend  = rend;
                tmpallinfo.mallele  = "";
                molAlignMarkerAllele.insert(std::pair<unsigned long, READINFO>(rstart, tmpallinfo));
                //
                rsp_itr ++;
                geti ++;
            }
            //
            map<unsigned long, string> alleleinfo;     // <pos_sorted, from>: only in given marker list
            map<unsigned long, double> allelefreq;     // <pos_sorted, af>
            map<unsigned long, string> allelecov;      // <pos_sorted, cov>
            // find allele info at markers of current molecule with the above range
            string bckey = lineinfo[0];
            std::pair <multimap<string, iMARKER>::iterator, 
                       multimap<string, iMARKER>::iterator> mkritr_region;
            mkritr_region = moleculeSet.equal_range(bckey);
            if(mkritr_region.first == mkritr_region.second)
            {   
                // current molecule does not cover any markers - skip.
                continue;
            }
            multimap<string, iMARKER>::iterator mkritr;
            for (mkritr=mkritr_region.first; mkritr!=mkritr_region.second; ++mkritr)
            {            
                iMARKER mkrtmp = (*mkritr).second;
                // firstly check if marker position covered by any alignments 2018-05-23
                bool rcmchecking = false;
                if(rcm==true && molReadCoveredPos.find(mkrtmp.pos) != molReadCoveredPos.end())
                {
                    rcmchecking = true;
                }
                if(rcm==false || rcmchecking==true)
                {    // record info at given markers
                    if(mkrtmp.pos>=msta && mkrtmp.pos<=mend && mkrtmp.isGiven==true) 
                    {
                        if(alleleinfo.find(mkrtmp.pos) != alleleinfo.end())
                        {
                            if(alleleinfo[mkrtmp.pos].compare(mkrtmp.from) != 0)
                            {
                                // this should not happen as het-cov has been checked above.
                                cout << "   Error: marker at " << bckey << ": " << mkrtmp.pos
                                     << " covering het but not filtered out. I need to check 1!! " << endl;
                                return false;
                            }
                        }    
                        else
                        {
                            alleleinfo.insert(std::pair<unsigned long, string>(mkrtmp.pos, mkrtmp.from));
                            allelefreq.insert(std::pair<unsigned long, double>(mkrtmp.pos, mkrtmp.af));
                            allelecov.insert(std::pair<unsigned long, string>(mkrtmp.pos, mkrtmp.cov));
                        }
                    }
                }
            }
            assert(allelefreq.size() == alleleinfo.size());
            assert(allelecov.size()  == alleleinfo.size());
            if(alleleinfo.size() == 0)
            {
                // no marker on this molecule
                continue;
            }
            
            //
            total_check ++;
            
            // build artificial sequence for current molecule
            unsigned long firstpos = 0; // first marker pos
            unsigned long lastpos  = 0; // last  marker pos
            bool firstposb = true;            
            string artSeq("");
            std::stringstream sspos;
            sspos.str("");
            map<unsigned long, string>::iterator alitr;
            map<unsigned long, string>::iterator alitr_end;
            alitr     = alleleinfo.begin();
            alitr_end = alleleinfo.end();
            map<unsigned long, string> toRemoveMarker; // those within the molecule but not covered by any reads
            while(alitr != alitr_end)
            {
                bool currentMmkrAdded = false;
                //
                std::stringstream sstmp;
                sstmp.str("");
                sstmp << (*alitr).first;
                //
                bool markerCoveredbyRead = false;
                // new 2018-05-23: reads covering same markers can happen
                multimap<unsigned long, READINFO>::iterator xitr;
                multimap<unsigned long, READINFO>::iterator xitr_end;
                xitr     = molAlignMarkerAllele.begin();
                xitr_end = molAlignMarkerAllele.end();
                while(xitr != xitr_end)
                {
                    unsigned long spansta = (*xitr).first;
                    unsigned long spanend = (*xitr).second.spanend;
                    
                    if(spansta<=(*alitr).first && (*alitr).first<=spanend)
                    {
                        markerCoveredbyRead = true;
                        if(currentMmkrAdded == false)
                        {
                            artSeq   += (*alitr).second;
                            currentMmkrAdded = true; // avoid adding mut allele multiple times
                        }
                        //                
                        if((*xitr).second.mallele.size()>0)
                        {
                            //string tmpstr(",\t----------\n         ");
                            string tmpstr(",");
                            (*xitr).second.mallele += (tmpstr + sstmp.str() + "\t" + (*alitr).second);
                        }
                        else
                        {
                            (*xitr).second.mallele += (         sstmp.str() + "\t" + (*alitr).second);
                        } 
                        if(firstposb)
                        {
                            firstpos = (*alitr).first;
                            lastpos  = (*alitr).first;
                            firstposb= false;
                        }     
                        else
                        {               
                            if((*alitr).first < firstpos)
                            {
                                firstpos = (*alitr).first;
                            }
                            else
                            if((*alitr).first > lastpos)
                            {
                                lastpos  = (*alitr).first;
                            }
                            else ;
                        }                                             
                    }
                    xitr ++;
                }
                //
                if(markerCoveredbyRead == false)
                {
                    // not covered by any read on this molecule: marker info needs to be removed
                    toRemoveMarker.insert(std::pair<unsigned long, string>((*alitr).first, "NF"));
                }
                //
                alitr ++;
            }
            // remove marker info that are not covered by any reads on the current molecule
            map<unsigned long, string>::iterator cl_itr;
            map<unsigned long, string>::iterator cl_itr_end;
            cl_itr     = toRemoveMarker.begin();
            cl_itr_end = toRemoveMarker.end();
            while(cl_itr != cl_itr_end)
            {
                unsigned long clkey = (*cl_itr).first;
                alleleinfo.erase(clkey);
                allelefreq.erase(clkey);
                allelecov.erase(clkey);
                cl_itr ++;
            }
            //
            multimap<unsigned long, READINFO>::iterator allitr;
            multimap<unsigned long, READINFO>::iterator allitr_end;
            allitr     = molAlignMarkerAllele.begin();
            allitr_end = molAlignMarkerAllele.end();
            while(allitr != allitr_end)
            {
                READINFO tmpallinfo = (*allitr).second;
                if(tmpallinfo.mallele.size()==0)
                {
                    tmpallinfo.mallele = "00000000\t0";
                }
                vector<string> mkrinfo = split_string(tmpallinfo.mallele, ',');
                vector<string>::iterator mitr     = mkrinfo.begin();
                vector<string>::iterator mitr_end = mkrinfo.end();
                while(mitr != mitr_end)
                {
                    sspos << "         " << (*mitr)             << "\t" 
                                         << tmpallinfo.spansta  << "\t" << tmpallinfo.spanend << "\t"
                                         << tmpallinfo.mateid   << "\t"
                                         << tmpallinfo.readname << "\n";
                    mitr ++;
                }
                allitr ++;
            }
            
            // check if there is breakpoint for current molecule using artificial sequence
            int allele1 = std::count(artSeq.begin(), artSeq.end(), '1');
            int allele2 = std::count(artSeq.begin(), artSeq.end(), '2');
            if(allele1>=min_marker && allele2>=min_marker )
            {
                cout << "   check: bckey " << bckey << " for current molecule " << msta << " to " <<  mend << ":" << endl;
                cout << "        : firstmkrpos->lastmkrpos = " << firstpos << "->" << lastpos << endl;
                cout << "        : artificial seq for current molecule " << artSeq << endl;
                cout << sspos.str() << endl;
                if(readwarning)
                {
                    cout << "      Warning: there were reads with same start position, last kept only. " << endl;
                }
                
                unsigned long leftp; // real
                unsigned long rightp;// real
                int           bpPos; // position breaking artSeq
                double        score;        
                if(get_break_pos_of_artSeq(artSeq, alleleinfo, min_marker, &leftp, &rightp, &bpPos, &score))
                {
                    // check size of potential CO interval
                    if(rightp - leftp + 1 <= maxCOSize)
                    {
                        unsigned long leftspan  = leftp   - firstpos + 1;
                        unsigned long rightspan = lastpos - rightp   + 1;
                        if(leftspan>=min_span  && rightspan>=min_span ) // caution - fp and fn    
                        {
                            if(score > min_score)
                            {
                                /* check read density within potential CO interval:
                                   if the size of the CO interval is too large, we still expect
                                      a number of reads within this interval, although
                                      they may not cover any markers.
                                   if the size of the CO interval is too small, however, it can be 
                                      due to the false mergeing of molecules.
                                   
                                   500 << bp_per_read_in_co <= 5000 bp - can be defined as good candidates in our settings.
                                   
                                   Note: this should be set according to real sequencing coverage on molecules.
                                */
                               int  bp_per_read_in_co         = -1;      // read density within potential CO interval
                               bool bp_boundary_same_fragment = false;   // whether reads covering co start/end from the same fragment
                               int  min_RR_dist               = 9999999; // minimum R1-R1 or R2-R2 distance within the molecule
                               int  max_RR_dist               = 0;       // maximum ...
                               int  n_readLeft                = 0;       // number of uniq reads ids on CO left (including non-marker covering)
                               int  n_readRight               = 0;       // number of uniq reads ids on CO right (including non-marker covering)
                               int  max_CORR_dist             = 0;       // maximum R1-R1 distance within CO interval
                                
                               if(!check_co_read_density(leftp,
                                                         rightp,
                                                         lineinfo[8], 
                                                         lineinfo[9],
                                                         lineinfo[10],
                                                         &bp_per_read_in_co,
                                                         &bp_boundary_same_fragment,
                                                         &min_RR_dist,
                                                         &max_RR_dist,
                                                         &n_readLeft,
                                                         &n_readRight,
                                                         &max_CORR_dist))
                               {
                                   return false;
                               }
                               // find allele freq and cov at markers
                               std::stringstream artSeqCov;
                               artSeqCov.str("");            
                               std::stringstream artSeqAF;
                               artSeqAF.str("");
                               map<unsigned long, string>::iterator alitri;
                               map<unsigned long, string>::iterator alitri_end;
                               alitri     = alleleinfo.begin();
                               alitri_end = alleleinfo.end();
                               int ibp    = 0;
                               while(alitri != alitri_end)
                               {
                                   //
                                   if(artSeqAF.str().size()>0 && ibp==bpPos)
                                   {
                                       artSeqAF << ";";
                                   }
                                   else
                                   if(artSeqAF.str().size()>0)
                                   {
                                       artSeqAF << ",";
                                   }
                                   assert(allelefreq.find((*alitri).first) != allelefreq.end());
                                   artSeqAF << allelefreq[(*alitri).first];
                                   //
                                   if(artSeqCov.str().size()>0 && ibp==bpPos)
                                   {
                                       artSeqCov << ";";
                                   }
                                   else
                                   if(artSeqCov.str().size()>0)
                                   {
                                       artSeqCov << ",";
                                   }
                                   assert(allelecov.find((*alitri).first) != allelecov.end());
                                   artSeqCov << allelecov[(*alitri).first];                
                                   //
                                   ibp ++;
                                   alitri ++;
                               }
                               // set minCORDist as 600 
                               if( (minCORDist<=bp_per_read_in_co && bp_per_read_in_co<=maxCORDist) ||
                                   (minCORDist>=bp_per_read_in_co && n_readLeft>=minRead && n_readRight>=minRead)
                                 )
                               {
                                    if( 1 )
                                    {
                                        //vector<string> chrbarcode = split_string(lineinfo[0], '#');
                                        std::stringstream sstmp;
                                        sstmp.str("");
                                        sstmp << chrbarcode[0]             << "\t"    
                                              << leftp                     << "\t"
                                              << rightp                    << "\t"
                                              << score                     << "\t"
                                              << artSeq.substr(0, bpPos)   << "\t" 
                                              << artSeq.substr(bpPos)      << "\t"
                                              << chrbarcode[1]             << "\t"
                                              << lastpos - firstpos + 1    << "\t"
                                              << leftp   - firstpos + 1    << "\t"
                                              << lastpos - rightp   + 1    << "\t"
                                              << bp_per_read_in_co         << "\t"
                                              << bp_boundary_same_fragment << "\t"
                                              << min_RR_dist               << "\t"
                                              << max_RR_dist               << "\t"
                                              << molelen                   << "\t"
                                              << readnum                   << "\t"
                                              << molecov                   << "\t"
                                              << artSeqAF.str()            << "\t"
                                              << artSeqCov.str()           << "\t"
                                              << rightp - leftp + 1        << "\t"
                                              << n_readLeft                << "\t"
                                              << n_readRight               << "\t"
                                              << max_CORR_dist;
                                        rawfp <<  sstmp.str() << endl;
                                        
                                        // caution filter out larger CO intervals
                                        // if(rightp - leftp < 15000)
                                        string aord = check_allele_order(artSeq.substr(0, bpPos), artSeq.substr(bpPos));
                                        insert_raw_bp(&rawbpmap, 
                                                      chrbarcode[0], 
                                                      leftp, 
                                                      rightp, 
                                                      score,
                                                      aord, 
                                                      sstmp.str());  
                                        cout << "        bp: " << leftp << " -> " << rightp 
                                             << "\t read density in interval = "  << bp_per_read_in_co << " bp/read" << endl;
                                        cout << "        allele frequency at markers: " << artSeqAF.str()  << endl;
                                        cout << "        allele coverage  at markers: " << artSeqCov.str() << endl;
                                        // collect detailed molecule info with a CO prediction
                                        std::stringstream coinfo;
                                        coinfo.str("");
                                        cout << "   Info: prepare file for " << chrbarcode[0]+"#" << leftp << "#" << rightp << endl;         
                                        coinfo << outPrefix+ipre.str() << "/" 
                                               << chrbarcode[0]        << "_" 
                                               << leftp                << "_" 
                                               << rightp               << "_" 
                                               << aord                 << "_"
                                               << chrbarcode[1]        << ".txt\0";
                                        cout << "   Info: details of molecule recorded in file " << coinfo.str() << endl;                                           
                                        ofstream tmpofp;
                                        tmpofp.open((coinfo.str()).c_str(), ios::out | ios::app);
                                        if(!tmpofp.good())
                                        {
                                            cout   << "   Error: cannot open file " << coinfo.str() << endl;
                                            return false;
                                        }
                                        else
                                        {
                                            tmpofp << "# " << sstmp.str() << endl; 
                                            tmpofp << "# bckey " << bckey << " for current molecule " << msta << " to " <<  mend << ":" << endl;
                                            tmpofp << sspos.str() << endl;
                                            tmpofp.close();
                                        }
                                    }
                                    else
                                    {
                                        cout << "        bp: " << leftp << " -> " << rightp 
                                             << "\tsame fragment but not sufficient uniq reads on each side " 
                                             << endl;
                                    }
                                }
                                else
                                {
                                    cout << "        bp: " << leftp << " -> " << rightp 
                                         << "\t read density in interval = " << bp_per_read_in_co << " bp/read "
                                         << " not in given range [" << minCORDist << ", " << maxCORDist << "]" << endl;
                                    cout << "        allele frequency at markers: " << artSeqAF.str() << endl;
                                    cout << "        allele coverage  at markers: " << artSeqCov.str() << endl;
                                }                      
                            }
                            else
                            {
                                cout << "        bp: " << leftp << " -> " << rightp << "\tlow score=" << score << endl;                
                            }
                        }
                        else
                        {
                            cout << "        bp: " << leftp << " -> " << rightp << "\tsmall span" << endl;         
                        }
                        
                    } // end of CO size checking
                    else
                    {
                        cout << "        larger interval size: " << rightp - leftp + 1 << " bp." << endl;
                    }
                } // end of if(get_break_pos_of_artSeq(artSeq, alleleinfo, min_marker, &leftp, &rightp, &bpPos, &score))
                else
                {
                    cout << "   Error: failed in scoring molecule sequence. " << endl;
                    return false;
                }
            } 
            else
            {
                // length of Col-/Ler-molecules
                if(outLerCol && allele1*1.0/(allele1+allele2) >= 0.99) // Col-molecule
                {
                    colMolefp << bckey << "\t" << lineinfo[1] << "\t" << lineinfo[3] << "\t" << artSeq << endl;
                }
                else
                if(outLerCol && allele2*1.0/(allele1+allele2) >= 0.99) // Ler-molecule                
                {
                    lerMolefp << bckey << "\t" << lineinfo[1] << "\t" << lineinfo[3] << "\t" << artSeq << endl;
                }
                else
                {
                    ;
                }
            }
        }
        // close file of raw molecules for chr
        molfp.close();
        
        vchritr ++;
    }
    // close files
    goodmolefp.close();
    if(outLerCol)
    {
        colMolefp.close();
        lerMolefp.close();
    }    
    cout << "   Info: total number of >=1 kb molecules checked: " 
         << numbTotalMol     
         << ", with "
         << numbTotalMolRead 
         << " reads. " 
         << endl;
    cout << "   Warning: " 
         << numbad        
         << " \"false\" molecules (based on RN, ML and HET), where " 
         << endl;    
    cout << "               " 
         << numbadRN      
         << " did not pass checking on RN <=" 
         << maxReadNum 
         << "; "    
         << endl;
    cout << "               " 
         << numbadML      
         << " did not pass checking on ML >=" 
         << 2*min_span 
         << " bp & <=" 
         << maxMoleLen 
         << " bp;"  
         << endl;
    cout << "               " 
         << numbadHET     
         << " did not pass checking on HET (among those passing checkings on RN and ML); "
         << endl;
    cout << "               and in total there were "
         << badBarcodeCov 
         << " barcodes cover both mutant and reference alleles (before RN and ML checking). " 
         << endl;
    cout << "   Info: total number of \"real\" molecules passed checking on RN, ML and HET: " 
         << numbGoodTotalMol 
         << ", with " 
         << numbGoodTotalMolRead 
         << " (" << setprecision(2)  << fixed
         << (float)numbGoodTotalMolRead/numbTotalMolRead *100
         << "%) reads. " 
         << endl;
    // output chr:left-coordinate-sorted raw break points: 
    // left-coordinate sorted raw break points
    string sortedrawbpfile = outPrefix + ipre.str() +"_raw_BP_sorted.txt";
    if(!output_co(rawbpmap, sortedrawbpfile))
    {
        return false;
    }
    // post-process predicated break points: final bp file
    string finalbpfile = outPrefix + ipre.str() +"_raw_BP_final.txt";
    if(!postprocess_bpPrediction(rawbpmap, finalbpfile))
    {
        cout << "   Error: processing raw break points failed. " << endl;
        return false;
    }
    if( !output_real_molecule_len_distribution(molecule_len_control, outPrefix+ipre.str()) )
    {
        cout << "   Error: output length distribution of real molecules failed. " << endl;
        return false;
    }
    // 
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    return true;
}
// sub
bool output_real_molecule_len_distribution(map<int, unsigned long> molecule_len_control, string outPrefix)
{
    // all "real" molecule passing checkings on ML, RN and HET
    // among which COs were reported.
    string outMLenFile = outPrefix + "_realMoLen_stat.txt";    
    ofstream ofp;
    ofp.open(outMLenFile.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open " << outMLenFile << " to write data. " << endl;
        return false;
    }
    map<int, unsigned long>::iterator mlitr;
    map<int, unsigned long>::iterator mlitr_end;
    mlitr     = molecule_len_control.begin();
    mlitr_end = molecule_len_control.end();
    ofp << "#molecule_len_kb\tmolecule_num_in_bin" << endl;
    while(mlitr != mlitr_end)
    {
        ofp << (*mlitr).first  << "\t" 
            << (*mlitr).second << endl;
        mlitr ++;
    }
    ofp.close(); 
    return true;
}
//    
bool check_mole_local_base_cov(unsigned long   mole_start,
                               unsigned long   mole_end,
                               string          read_align_info,
                               unsigned long   win_size,
                               unsigned long   win_step,
                               double*         min_win_base_cov,
                               double*         max_win_base_cov,
                               double*         mean_win_base_cov,
                               double*         sd_win_base_cov)
{
    
    return true;
} 
//
bool build_artSeq(map<unsigned long, string> alleleinfo, 
                  unsigned long regionBegin, 
                  unsigned long regionEnd, 
                  string*       artSeq,
                  string*       artSeqCheckingStr)
{
    // build artificial allele sequence within a given region at all possible markers.
    if(alleleinfo.size()==0) return false;
    std::stringstream posAllele;
    posAllele.str("");
    (*artSeq).clear();
    (*artSeqCheckingStr).clear();
    map<unsigned long, string>::iterator alitr;
    map<unsigned long, string>::iterator alitr_end;
    alitr     = alleleinfo.begin();
    alitr_end = alleleinfo.end();
    while(alitr != alitr_end)
    {
        if( (*alitr).first>=regionBegin && (*alitr).first<=regionEnd)
        {
            (*artSeq)   += (*alitr).second;
            posAllele << "        " << (*alitr).first << ":" << (*alitr).second << "\n";
        }
        alitr ++;
    }
    (*artSeqCheckingStr) = posAllele.str();
    return true;
}
//
bool check_co_read_density(unsigned long       co_start,
                           unsigned long       co_end,
                           string              read_align_info,
                           string              read_order_info,
                           string              read_id_info,
                           int*                bp_per_read_in_co,
                           bool*               bp_boundary_same_fragment,
                           int*                min_RR_dist,
                           int*                max_RR_dist,
                           int*                n_readLeft,
                           int*                n_readRight,
                           int*                max_CORR_dist)
{
    /* this function checks read density in a potential CO interval
       how many bp per read within this interval
       -The larger/"smaller" the value, the higher chance that 
        this molecule is falsely merged from different molecules.
       
       Example molecule with info:
       
	1#ACGCGTGGTCCAGTAT  8229787 8277108 47422   0.08    30  R   8277208
 
 	# span of reads: read_align_info
	8229787,8229914;8229933,8230022;8234823,8234973;8235084,8235141;
	8235256,8235406;8235546,8235673;8237785,8237913;8238134,8238284;
	8245266,8245416;8245555,8245682;8250091,8250218;8250471,8250611;
	8252479,8252629;8252479,8252629;8252711,8252807;8252714,8252807;
	8256708,8256858;8257166,8257293;8262644,8262803;8262929,8263056;
	8264841,8264991;8265091,8265218;8266796,8266946;8266809,8266864;
	8267521,8267644;8267724,8267873;8272690,8272817;8273056,8273206;
	8276756,8276883;8277108,8277208 

	# order in read pairs: read_order_info
	1,2,2,1,2,1,1,2,2,1,1,2,2,2,1,1,2,1,2,1,2,1,2,1,1,2,1,2,1,2 

	# id of reads: read_id_info
	J00137:105:HNFWJBBXX:3:1213:4685:18616,J00137:105:HNFWJBBXX:3:1213:4685:18616,
	J00137:102:HMYVWBBXX:7:2102:6482:8400,J00137:102:HMYVWBBXX:7:2102:6482:8400,
	J00137:105:HNFWJBBXX:3:2225:2950:3266,J00137:105:HNFWJBBXX:3:2225:2950:3266,
	J00137:108:HNGC7BBXX:3:1211:16072:22133,J00137:108:HNGC7BBXX:3:1211:16072:22133,
	J00137:105:HNFWJBBXX:3:1109:7314:28973,J00137:105:HNFWJBBXX:3:1109:7314:28973,
	J00137:108:HNGC7BBXX:3:2103:14874:10247,J00137:108:HNGC7BBXX:3:2103:14874:10247,
	J00137:102:HMYVWBBXX:7:2127:1692:16559,J00137:102:HMYVWBBXX:7:2127:1895:16665,
	J00137:102:HMYVWBBXX:7:2127:1692:16559,J00137:102:HMYVWBBXX:7:2127:1895:16665,
	J00137:105:HNFWJBBXX:3:1114:28107:1525,J00137:105:HNFWJBBXX:3:1114:28107:1525,
	J00137:102:HMYVWBBXX:8:2124:19471:30274,J00137:102:HMYVWBBXX:8:2124:19471:30274,
	J00137:108:HNGC7BBXX:3:1210:19705:2897,J00137:108:HNGC7BBXX:3:1210:19705:2897,
	J00137:102:HMYVWBBXX:7:2106:30655:43251,J00137:102:HMYVWBBXX:7:2106:30655:43251,
	J00137:108:HNGC7BBXX:3:1212:11363:38029,J00137:108:HNGC7BBXX:3:1212:11363:38029,
	J00137:105:HNFWJBBXX:3:2204:5284:22854,J00137:105:HNFWJBBXX:3:2204:5284:22854,
	J00137:105:HNFWJBBXX:3:1104:14143:15451,J00137:105:HNFWJBBXX:3:1104:14143:15451
	
	bp predicted as: 8252716 -> 8262697, where CO start/end in reads at 8252714, 8262644
    */
    //
    if(co_end < co_start)
    {
        return false;
    }
    // return values
    *bp_per_read_in_co         = -1;     // read density within CO interval
    *bp_boundary_same_fragment = false;  // whether the read covering co start/end from the same fragment
    *min_RR_dist               = 999999; // minimum R1-R1 or R2-R2 distance within the molecule
    *max_RR_dist               = 0;      // maximum ...
    // task 1. check read density within the potential CO interval
    vector<string> respan  = split_string(read_align_info, ';'); // ; with ,    
    vector<string> reorder = split_string(read_order_info, ','); // ,
    vector<string> reid    = split_string(read_id_info,    ','); // ,
    assert( respan.size() == reorder.size() );
    assert( respan.size() == reid.size() ); 
    vector<string>::iterator ritr;
    vector<string>::iterator ritr_end;
    ritr     = respan.begin();
    ritr_end = respan.end();
    int n_read = 0;                    // number of reads in interval
    int i_read = -1;                   // to satisfy index start from 0 in vector
    unsigned long sta_aligned_pos = 0; // last  read covering start of interval
    unsigned long end_aligned_pos = 0; // first read covering end   of interval
    string        sta_read_name("");
    string        end_read_name("");
    string        sta_read_order("");
    string        end_read_order("");
    //
    map<string, int> readLeft;
    map<string, int> readRight;
    //
    unsigned long lastReadStaInCO     = 0; // <read_align_sta, read_align_end>
    unsigned long maxIntervalSizeInCo = 0;
    //
    while(ritr != ritr_end)
    {
        vector<string> ralign = split_string(*ritr, ',');
        
        unsigned long rsta = strtoul(ralign[0].c_str(), NULL, 0);
        unsigned long rend = strtoul(ralign[1].c_str(), NULL, 0);
        
        i_read ++; // first valid: i_read = 0
        
        // skip checking reads in front of the interval
        if(rend < co_start)
        {
            if(readLeft.find(reid[i_read]) == readLeft.end())
            {
                readLeft.insert(std::pair<string, int>(reid[i_read], 1));
            }
            ritr ++;
            continue;
        }
        //
        if(rsta<=co_start && co_start<=rend)
        {
            sta_aligned_pos = rsta; // take the last  read statisfying this
            sta_read_order  = reorder[i_read];            
            sta_read_name   = reid[i_read];
        }
        //
        if(end_aligned_pos==0 && rsta<=co_end && co_end<=rend)
        {
            end_aligned_pos = rsta; // take the first read satsifying this
            end_read_order  = reorder[i_read];            
            end_read_name   = reid[i_read];        
        }
        // number of reads within interval
        if(co_start<rsta && rend<co_end)
        {
            n_read ++;
        }
        // reads on/in CO interval
        if( (co_start<=rsta && rsta<=co_end) || (co_start<=rend && rend<=co_end) )
        {
            if(lastReadStaInCO != 0)
            {
                assert(rsta>=lastReadStaInCO);
                unsigned long this_dist = rsta - lastReadStaInCO;
                if(this_dist > maxIntervalSizeInCo) 
                {
                    maxIntervalSizeInCo = this_dist;
                }
            }
            lastReadStaInCO     = rsta;
        }        
        // reads after the interval
        if(rsta > co_end)
        {
            if(readRight.find(reid[i_read]) == readRight.end())
            {
                readRight.insert(std::pair<string, int>(reid[i_read], 1));
            }
        }
        ritr ++;
    }
    //
    *max_CORR_dist = maxIntervalSizeInCo;
    // exclude read ids at start and end of CO
    if(readLeft.find(sta_read_name) != readLeft.end())
    {
        readLeft.erase(sta_read_name);
    }
    *n_readLeft = readLeft.size();
    if(readRight.find(end_read_name) != readRight.end())
    {
        readRight.erase(end_read_name);
    }
    *n_readRight = readRight.size();
    // e.g., one read r2 within interval: co_sta_in_r1----r2------co_end_in_r3, that is, n_read = 1,
    //      => density = (co_end - co_sta) / (n_read+1)
    n_read += 1; 
    *bp_per_read_in_co = (co_end - co_start)/n_read; // check::TODO: use end_last_read - sta_first_read :: full span
    //
    string frag(" different ");
    if(end_read_name.compare(sta_read_name) == 0)
    {
        *bp_boundary_same_fragment = true;
        frag                       = (string)" same ";
    }
    //
    if(sta_read_name.size() != 0 && end_read_name.size() != 0)
    {
        cout << "         CO starts/ends in reads at " 
             << sta_aligned_pos   << ", " << end_aligned_pos << "\n"
             << "          with " << frag << " read ids: " 
             << sta_read_name     << ", " << end_read_name
             << endl;
    }
    else // reads in molecule missing compared to vcf, possibly due to secondary alignment etc.
    if(sta_read_name.size() != 0 && end_read_name.size() == 0)
    {
        cout << "         CO starts in read at " 
             << sta_aligned_pos   << ", CO-end read discarded during bam conversion. " << "\n"
             << "          with " << frag << " read id: " 
             << sta_read_name     << ","
             << endl;    
    }
    else
    if(sta_read_name.size() == 0 && end_read_name.size() != 0)
    {
        cout << "         CO ends in read at " 
             << "CO-start read discarded during bam conversion" << ", " << end_aligned_pos << "\n"
             << "          with " << frag << " read id: " 
             << ", "              << end_read_name
             << endl;    
    }
    else
    {
        cout << "         CO starts/ends in reads -- both discarded during bam conversion. " 
             << endl;
    }

    // task 2. check read density within the potential-CO related molecule
    unsigned long RR_dist_total  = 0;
    int           RR_dist_num    = 0;
    unsigned long lastR1AlignSta = 0; // for calculating R1 distance
    unsigned long lastR1AlignEnd = 0; // ..
    unsigned long lastR2AlignSta = 0; // for calculating R2 distance
    unsigned long lastR2AlignEnd = 0; // ..  
    ritr        = respan.begin();
    ritr_end    = respan.end();
    i_read      = -1;
    while(ritr != ritr_end)
    {
        vector<string> ralign = split_string(*ritr, ',');
        unsigned long rsta = strtoul(ralign[0].c_str(), NULL, 0);
        unsigned long rend = strtoul(ralign[1].c_str(), NULL, 0);
        // the i_read-th read: 0,1,...,n-1
        i_read ++;
        string iorder = reorder[i_read];
        
        if(i_read == 0)
        {
            if(iorder.compare("1") == 0)
            {
                lastR1AlignSta = rsta;
                lastR1AlignEnd = rend;
            }
            else
            {
                lastR2AlignSta = rsta;
                lastR2AlignEnd = rend;                
            }
        }
        else
        {
            if(iorder.compare("1") == 0) // current is R1
            {
                if(lastR1AlignSta > 0)   // there is last R1
                {
                    // the dist where a read takes within a molecule: always >= 0
                    int i_dist   = (rsta>lastR1AlignSta)?(rsta-lastR1AlignSta):((lastR1AlignSta-rsta)*(-1));
                    *min_RR_dist = (*min_RR_dist>i_dist)?i_dist:*min_RR_dist;
                    *max_RR_dist = (*max_RR_dist<i_dist)?i_dist:*max_RR_dist;
                    RR_dist_total+= i_dist>0?i_dist:(-1)*i_dist; // absolute dist total
                    RR_dist_num ++;
                }
                lastR1AlignSta  = rsta;
                lastR1AlignEnd  = rend;
            }
            else                         // current is R2
            {
                if(lastR2AlignSta > 0)   // there is last R2
                {
                    // minus exists: meaning overlapped RR reads
                    int i_dist  = (rsta>lastR2AlignSta)?(rsta-lastR2AlignSta):((lastR2AlignSta-rsta)*(-1));
                    *min_RR_dist = (*min_RR_dist>i_dist)?i_dist:*min_RR_dist;
                    *max_RR_dist = (*max_RR_dist<i_dist)?i_dist:*max_RR_dist;
                    RR_dist_total+= i_dist>0?i_dist:(-1)*i_dist; // absolute dist total
                    RR_dist_num ++;              
                }                
                lastR2AlignSta  = rsta;
                lastR2AlignEnd  = rend;                
            }          
        }
        
        ritr ++;
    }

    return true;
}
//
bool check_candidate_cosite(multimap<string, std::pair<unsigned long, unsigned long> >* badMolecule,
                            string mkey,
                            unsigned long mstart,
                            unsigned long mend)
{
    bool passchecking = true; // if not in bad molecule set, passing checking
    std::pair <multimap<string, std::pair<unsigned long, unsigned long> >::iterator, 
               multimap<string, std::pair<unsigned long, unsigned long> >::iterator> clitr_region;
    clitr_region = (*badMolecule).equal_range(mkey);
    multimap<string, std::pair<unsigned long, unsigned long> >::iterator clitr;
    for (clitr=clitr_region.first; clitr!=clitr_region.second; ++clitr)
    {
        std::pair<unsigned long, unsigned long> span = (*clitr).second;
        if( (mstart>=span.first && mstart<=span.second) ||
            (mend>=span.first   && mend<=span.second) )
        {
            // caution: might be too stringent here: check overlapping size then?
            passchecking = false;
            break;
        }
    }
    return passchecking; 
}
//
bool get_bad_molecule(const char* file, 
                      int         maxReadNum,
                      int         maxMoleLen,
                      multimap<string, std::pair<unsigned long, unsigned long> >* badMolecule)
{
    // 1. if molecule length > cutoff,      collect it as bad
    // 2. if molecule read number > cutoff, collect it as bad
    // 3. if rNumber > mLen/mLen_cutoff * rNumber_cutoff, collect it as bad - 2018-07-30
    // molecule table: statistics of all moleclues: generated by 'DrLink molecule'
    // format: chr#barcode    first_aligned   last_aligned    molecule_len    molecule_cov    read_num    Uni_flag
    igzstream fp;
    fp.open(file);
    if(!fp.good())
    {
        cout << "   Error: cannot open read number distribution " << file << endl;
        return false;
    }
    else
    {
        cout << "   Info: collecting bad molecules (with larger read number or size)... " << endl;
    }
    unsigned long numibad = 0;
    while(fp.good())
    {
        string line("");
        getline(fp, line);
        if(line.size()==0 || line[0]=='#') continue;
        
        vector<string> lineinfo = split_string(line, '\t');
        if(line.size() < 7) continue;
        
        int readnum = atoi(lineinfo[5].c_str());
        int molelen = atoi(lineinfo[3].c_str()); 
        if( (readnum>maxReadNum || molelen>maxMoleLen) ||  
            (readnum*maxMoleLen > maxReadNum*molelen)
          )// potential overlapping molecules
        {
            std::pair <unsigned long, unsigned long> span;
            span.first  = strtoul(lineinfo[1].c_str(), NULL, 0);
            span.second = strtoul(lineinfo[2].c_str(), NULL, 0);
            // chr+"#"+barcode => (sta, end) : allow multi-aligns with the same barcode.
            (*badMolecule).insert(std::pair<string, std::pair<unsigned long, unsigned long> >(lineinfo[0], span));  
            numibad ++;          
        }
    }
    fp.close();
    cout << "   Info: number of collected bad molecules with read number >  " << maxReadNum 
         << " or molecule length > "                                          << maxMoleLen << " bp: "
         << numibad << endl;    
    return true;
}
//
bool get_moleStat_distribution(const char* file, 
                              double      maxPercent, 
                              int*        maxStat,
                              string      lenORreadnum)
{
    // this function finds out a cutoff for read_number/length of molecules for determining some of them as bad ones.
    // lenORreadnum is either "readNum" or "moleLen"
    cout << "   Info: checking " << lenORreadnum << " distribution ... " << endl;
    fstream fp;
    fp.open(file, ios::in);
    if(!fp.good())
    {
        cout << "   Error: cannot open read number distributio " << file << endl;
        return false;
    }
    vector<int>           xvalue;
    vector<unsigned long> yvalue;
    //
    int maxNumLen = 0; // if about length, in kb
    unsigned long maxcnt = 0;
    int xindex    = 0;
    int maxxindex = 0;
    char sepc = '\t';
    while(fp.good())
    {
        string line("");
        getline(fp, line);
        if(line.size()==0 || line[0]=='#') continue;
        
        if(sepc == '\t' && line.find("\t") == std::string::npos) 
        {
            sepc = ' ';
        }
        vector<string> lineinfo = split_string(line, sepc);
        if(line.size() < 2) continue;
        
        int           xxx = round(atof(lineinfo[0].c_str()));
        unsigned long yyy = strtoul(lineinfo[1].c_str(), NULL, 0);
        xindex ++;
        
        xvalue.push_back(xxx);
        yvalue.push_back(yyy);
        
        if(yyy > maxcnt)
        {
            maxcnt    = yyy;
            maxNumLen = xxx;
            maxxindex = xindex; // 1-based system; indeed ymax is at vector[maxxindex-1]
        }
    }
    // find a proper maximum statistic
    assert(maxxindex>=1);
    for(int i = maxxindex; i < xvalue.size(); i ++)
    {
        // caution: too stringent on molecule size
        if( (yvalue[i]*1.0)/(yvalue[maxxindex-1]*1.0)<maxPercent )
        {
            if(i == maxxindex)
            {
                // the one right next to max
                *maxStat = xvalue[i]; 
            }
            else
            {
                *maxStat = xvalue[i-1];
            }
            break;
        }
        else 
        if( lenORreadnum.compare("readNum")==0 && 
            (double)xvalue[i] > 1.5*(double)xvalue[maxxindex-1])
        {
            *maxStat = xvalue[i-1];
            break;
        }
        else 
        if( lenORreadnum.compare("moleLen")==0 && 
            (double)xvalue[i] > 1.5*(double)xvalue[maxxindex-1])
        {
            *maxStat = xvalue[i-1];
            break;
        }        
    }
    if(lenORreadnum.compare("readNum")==0)
    {
        *maxStat -= 1;
        cout << "   Info: maximum number of reads allowed for a good molecule: " 
             << *(maxStat)   << ";"    
             << endl;
        cout << "   Info: checking read number distribution done. "  
             << endl;
    }
    else
    if(lenORreadnum.compare("moleLen")==0)
    {
        *maxStat = (*maxStat) * 1000 + 1000; // into bp and complement 1kb: less stringent
        cout << "   Info: maximum length allowed for a good molecule: "          
             << *maxStat     
             << " bp;" 
             << endl;
        cout << "   Info: checking molecule length distribution done. "
             << endl;
    }   
    else
    {
        cout << "   Error: you should never reach here when checking readnum or length of molecules." 
             << endl;
        return false;
    } 
    return true;
}
//
bool get_break_pos_of_artSeq(string                     artSeq, 
                             map<unsigned long, string> artSeqPos, 
                             int                        min_marker,
                             unsigned long*             leftp, 
                             unsigned long*             rightp, 
                             int*                       bpPos,
                             double*                    score)
{
    *score  = 0.0;
    map<unsigned long, string>::iterator pitr;
    pitr = artSeqPos.begin();
    int i = 0;
    while(i < min_marker)
    {
        pitr ++;
        i ++;
    }
    for(size_t mi = min_marker; mi <= artSeq.size()-min_marker; mi ++)
    {
        string strLef = artSeq.substr(0, mi); // val @ pos of 0,...,mi-1
        string strRig = artSeq.substr(mi);    // val @ pos of mi,mi+1,...
        // current break point
        size_t c1Lef = std::count(strLef.begin(), strLef.end(), '1'); // L1
        size_t c2Lef = strLef.size()-c1Lef;                           // L2
        size_t c1Rig = std::count(strRig.begin(), strRig.end(), '1'); // R1
        size_t c2Rig = strRig.size()-c1Rig;                           // R2
        
        double af1Lef = (double)c1Lef/strLef.size();
        double af2Lef = (double)c2Lef/strLef.size();
        
        double af1Rig = (double)c1Rig/strRig.size();
        double af2Rig = (double)c2Rig/strRig.size();
        
        double score1 = af1Lef * af2Rig;
        double score2 = af2Lef * af1Rig;
        
        double scotmp = score1>score2?score1:score2;
        
        if(scotmp > *score)
        {
            *score  = scotmp;
            pitr --;
            *leftp  = (*pitr).first;
            pitr ++; 
            *rightp = (*pitr).first;
            *bpPos  = mi;
            if(*rightp == *leftp)
            {
                cout << "   artSeq  ="  << artSeq    << "\n"
                     << "   mi    ="  << *bpPos  << "\n"
                     << "   leftp ="  << *leftp  << "\n"
                     << "   rightp="  << *rightp << "\n" << endl;
                cout << "   positions related to artSeq: " << "\n   ";
                for(map<unsigned long, string>::iterator itr = artSeqPos.begin(); 
                    itr != artSeqPos.end(); 
                    itr ++) 
                {
                    cout << (*itr).first << " ";
                }
                cout << "\n   Error-must-fix: molecule separated at the same position -- impossible! " 
                     << endl;
                //return false;
            }
        }
        //
        pitr ++;
    }
    return true;
}
//
bool get_var_ref_barcode(string                formatInfo, 
                         string                barcodeString, 
                         vector<string>*       vbc, 
                         vector<string>*       rbc,
                         string                this_chr,
                         string                this_pos,
                         multimap<string, unsigned long>* HeteroBarcode)
{
    //    2017-10-30: filter interleaving reads
    //                check if there are same barcode-reads covering the same allele at some variants!!!
    //                more possibly, these reads are from molecules with the same chr#barcode but from different plants.
    
    // formatInfo   : GT:DP:RO:QR:AO:QA:GL:BX:PS - can be different but must contain 'BX' flag
    // barcodeString: real values corresponding to the above format
    // vbc          : barcodes related to variant   allele
    // rbc          : barcodes related to reference allele
    vector<string> fmt = split_string(formatInfo, ':');
    int pos = -1;
    for(int i = 0; i < fmt.size(); i ++)
    {
        if(fmt[i].compare("BX") == 0)
        {
            pos = i;
            break;
        }
    }
    if(pos < 0)
    {
        cout << "   Error: no BX found in FORMAT. " << endl;
        return false;
    }
    vector<string> bcode = split_string(barcodeString, ':');
    if(bcode.size() != fmt.size())
    {
        cout << "   Error: sizes of FORMAT and VALUE mismatched." << endl;
        return false;
    }
    // find barcodes related to all alleles
    string b_all     = bcode[pos];
    // find barcodes related to alternative allele, after ',': "bcr1-;bcr2-,bcv1-;bcv2-;bcv3-"
    size_t seRefAlt  = b_all.find(",");
    assert(seRefAlt != std::string::npos);    
    //
    string b_var     = b_all.substr(seRefAlt+1);
    if(b_var.size()>15)
    {
        (*vbc)           = split_string(b_var, ';');
        // remove redundant info
        vector<string>::iterator itr;
        vector<string>::iterator itr_end;
        itr        = (*vbc).begin();
        itr_end    = (*vbc).end();  
        while(itr != itr_end)
        {
            string tmp      = *itr;
            size_t trashpos = tmp.find("-");
            if(trashpos == std::string::npos)
            {
                trashpos = tmp.size();
            }
            if(trashpos > 0)
            (*itr)          = tmp.substr(0, trashpos);
            itr             ++;
        }
    }
    else
    {
        (*vbc).clear(); // possibly no variant barcode? no!
    }
    // find barcodes related to reference allele, before ',': "bcr1-;bcr2-,bcv1-;bcv2-;bcv3-"
    string b_ref     = b_all.substr(0, seRefAlt);
    if(b_ref.size()>15)
    {
        (*rbc)           = split_string(b_ref, ';');
        // remove redundant info
        vector<string>::iterator itr;
        vector<string>::iterator itr_end;        
        itr        = (*rbc).begin();
        itr_end    = (*rbc).end();
        while(itr != (*rbc).end()) // note that (*rbc).end() is changing, so itr_end cannot be used here!!!
        {
            string tmp      = *itr;
            size_t trashpos = tmp.find("-");
            if(trashpos == std::string::npos)
            {
                trashpos = tmp.size();
            }
            if(trashpos > 0)
            (*itr)          = tmp.substr(0, trashpos);
            vector<string>::iterator overlapitr = std::find((*vbc).begin(), (*vbc).end(), *itr);
            if(overlapitr != (*vbc).end())
            {
                if(badBarcodeCov == 0)
                {
                    cout << "   Info: barcode " <<  *itr << " of "     << barcodeString 
                         << " covers both reference and alternative. " << endl;
                    cout << "   Warning: molecules with such barcodes removed (a total number given in the end). " 
                         << endl;
                }
                else
                {
                    if(false)
                    {
                        cout << "   Info: barcode " <<  *itr << " of " 
                             << barcodeString << " covers both reference and alternative. " << endl;
                    }
                }
                // count
                badBarcodeCov ++;
                // this is the barcode for removing falsely merged molecule covering the heterozgous mkr pos
                //(*vbc).erase(overlapitr);
                //itr = (*rbc).erase(itr);
                unsigned long mkrposwithbadbc = strtoul(this_pos.c_str(), NULL, 0);
                (*HeteroBarcode).insert(std::pair<string, unsigned long>(this_chr+"#"+(*itr), mkrposwithbadbc));
            }
            // next barcode
            itr ++;
        }
    }
    else
    {
        (*rbc).clear(); // possibly no variant barcode? yes!
    }
    return true;
}
//
bool collect_marker_info(const char* file, string parent_id, map<string, string>* mkr)
{
    // type of markers: snps, 1-bp deletion and 1-bp insertion
    cout << "   Info: reading marker info from " << file << endl;
    fstream ifp;
    ifp.open(file);
    if(!ifp.good())
    {
        return false;
    }
    int num = 0;
    int raw = 0;
    int del = 0;
    int ins = 0;
    bool adjsnp = false; // SNPs directly next to each other
    int    nadj = 0;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        
        /* e.g.,
	        ks	1	164393	C	G  -- snp
		ks	1	164413	AC	TA -- twp snps
		ks	1	460940  -	A  -- insertion
		ks	1	468923  T	-  -- deletion
        */
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo[3].size() != lineinfo[4].size()) continue; // non-snps
        raw ++;
        
        if(lineinfo[3].size() > 1 && !adjsnp)
        {
            cout << "   Info: snp line will be separated, e.g, " << line << endl;
            cout << "       You will get a total number after reading all marker info. " << endl;
            adjsnp = true;
            nadj ++;
        }
        else
        if(lineinfo[3].size()>1)
        {
            nadj ++;
        }
        
        //  split "double/triple/..." snps
        for(int snpi = 0; snpi < lineinfo[3].size(); snpi ++)
        {
            string ref = lineinfo[3].substr(snpi, 1);
            string alt = lineinfo[4].substr(snpi, 1);
            //
            unsigned long pos = strtoul(lineinfo[2].c_str(), NULL, 0);
            pos += snpi;
            std::stringstream sspos;
            sspos.str("");
            sspos << pos;
            //
            string key = lineinfo[1] + "#" + sspos.str() + "#" + ref + "#" + alt; // "GA", "G-" or "-A"
            string val = parent_id;
            
            map<string, string>::iterator itr;
            map<string, string>::iterator itr_end;
            itr     = (*mkr).find(key);
            itr_end = (*mkr).end();
            
            if(itr == itr_end)
            {
                (*mkr).insert(std::pair<string, string>(key, val));
                num ++;
                if(alt=="-") del ++;
                if(ref=="-") ins ++;
            } 
            else
            {
                // cout << "    Warning: same marker for both parents: " << line << " (only recorded in first set)" << endl;
                if(parent_id.compare((*itr).second) != 0)
                {
                    (*mkr).erase(key); // removed
                }
            }
        }
    }
    ifp.close();
    cout << "   Info: recording " << num << " (of " << raw << ")" << " snp markers info done " << endl;  
    cout << "   Info: total number of \"neighboring\" snps: "     << nadj << endl;  
    cout << "   Info: total number of del markers: "              << del  << endl;
    cout << "   Info: total number of ins markers: "              << ins  << endl;
    if(parent_id.compare("2") == 0)
    {
        cout << "   Info: in total  " << (*mkr).size()  << " markers info recorded till now. " << endl; 
        cout << "   Info: markers appear in both marker sets have been removed from records. " << endl << endl;
    }
    return true;
}
//
bool read_chr_var(string vcffilename, map<string, string> knownMarker, multimap<string, iMARKER>* moleculeSet)
{
    // this function collects variants specific to chr from subfile_chr_x_longranger_variants.vcf.gz
    igzstream  fp;
    fp.open(vcffilename.c_str());
    if(!fp.good())
    {
        cout << "   Error: cannot open vcf file " 
             << vcffilename 
             << ". " 
             << endl;
        return false;
    }
    //
    unsigned long snpnum        = 0;
    unsigned long delnum        = 0;
    unsigned long insnum        = 0;
    unsigned long lineNum       = 0;
    unsigned long mkrNumChr     = 0;
    unsigned long snpnum_all    = 0;
    unsigned long delnum_all    = 0;
    unsigned long insnum_all    = 0;
    unsigned long mkrNumChr_all = 0;    
    cout << "   Info: only markers with coverage in [" 
         << minVarCov << "x, " << maxVarCov << "x] will be recorded (snps/1-bp indels considered). " 
         << endl;
    while(fp.good())
    {
        string line("");
        getline(fp, line);
        if(line.size()==0 || line[0]=='#' || line[0]=='*') continue;
               
        lineNum ++;
        if(lineNum%100000 == 0) cout << "   Info: " 
                                     << lineNum     << "th var of " 
                                     << vcffilename << "..." 
                                     << endl;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 10)
        {
            cout << "   Warning: skip line due to insufficient cell info: " << line << endl; 
            continue;
        }
        //
        if(lineinfo[4].find(",")!=std::string::npos)
        {
            continue; // multiple alternative allele; skip!
        }
        //
        if(lineinfo[6].compare("PASS") != 0)
        {
            // caution: USED: low-quality variation calls - 2018-06-07
            // check: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/vcf
            // continue;
        }
        // get barcodes related to marker; especially check het-barcodes before any filtering
        vector<string> varBarcodes;   
        vector<string> refBarcodes;         
        if(!get_var_ref_barcode(lineinfo[8], 
                                lineinfo[9], 
                                &varBarcodes, 
                                &refBarcodes, 
                                lineinfo[0], 
                                lineinfo[1], 
                                &HeteroBarcode))
        {
            cout << "   Warning: unexpected barcode info at " << line << endl;
            continue;
        }
        // key= chr#pos#ref#alt: snp: 1#1000#G#A, del: 1#1000#G#-, ins: 1#1000#-#G
        string key;
        string itype("snp");
        std::stringstream dss;
        dss.str("");        
        if(lineinfo[3].size()==2 && lineinfo[4].size()==1) // 1-bp del
        {
            unsigned long dpos = strtoul(lineinfo[1].c_str(), NULL, 0);
            dpos ++;
            dss << dpos;
            key   = lineinfo[0] + "#" + dss.str() + "#" + lineinfo[3].substr(1,1) + "#" + "-";
            itype = "del";
        }
        else
        if(lineinfo[3].size()==1 && lineinfo[4].size()==2) // 1-bp ins
        {
            dss << lineinfo[1];
            key   = lineinfo[0] + "#" + dss.str() + "#" + "-" + "#" + lineinfo[4].substr(1,1);
            itype = "ins";            
        }
        else
        if(lineinfo[3].size()==1 && lineinfo[4].size()==1) // snp
        {
            dss << lineinfo[1];
            key   = lineinfo[0] + "#" + dss.str() + "#" + lineinfo[3] + "#" + lineinfo[4];
        }
        else
        {
            // skip >=2bp non SNP variations/multiple alleles; // new 2018-06-16
            continue;
        }        
        //
        // all markers checked
        if(itype.compare("del")==0) delnum_all ++;
        else
        if(itype.compare("ins")==0) insnum_all ++;
        else
        if(itype.compare("snp")==0) snpnum_all ++;
        else ;
        mkrNumChr_all ++;        
        //
        bool isGiven = true;
        if(knownMarker.find(key) == knownMarker.end()) 
        {
            isGiven = false; // not in marker list; skip!
            continue;
        }        
        // check raw coverage on marker site - caution - only well-covered snps used as markers.
        size_t posDPsta = lineinfo[7].find("DP=");
        size_t posDPend = lineinfo[7].find(";", posDPsta);
        string posCov   = lineinfo[7].substr(posDPsta+3, posDPend - posDPsta - 3);
        
        // testing on allele coverage: 2018-04-11: sample A: avg 323+-323/4
        int allele_depth = strtoul(posCov.c_str(), NULL, 0);
        if(allele_depth<minVarCov || allele_depth>maxVarCov)
        {
            /*
            cout << "   Check: low-cov marker " << posCov << "x at " 
                 << lineinfo[0] << "\t" 
                 << lineinfo[1] << "\t" 
                 << lineinfo[3] << "\t"
                 << lineinfo[4] << endl;
            */
            continue;
        }
        // check observed allele frequency of markers: from AB <= SRF;SRR;SAF;SAR or RO/AO
        // SRF=55;SRR=62;SAF=43;SAR=53 ==> AB=(43+53)÷(43+53+55+62)
        double maf = 0; // allele frequency
        /*
        if(lineinfo[7].find("AB=") != std::string::npos)
        {
            size_t posABsta = lineinfo[7].find("AB=");
            size_t posABend = lineinfo[7].find(";", posABsta+1); 
            if(posABend == std::string::npos)
            {
                posABend = lineinfo[7].size()+1;
            }
            string ABstr    = lineinfo[7].substr(posABsta+3, posABend - posABsta - 3);
            maf             = atof(ABstr.c_str());           
        }
        else
        if(lineinfo[7].find("SRF=")!=std::string::npos &&
           lineinfo[7].find("SRR=")!=std::string::npos &&
           lineinfo[7].find("SAF=")!=std::string::npos &&
           lineinfo[7].find("SAR=")!=std::string::npos)
        {
            size_t posSRFsta = lineinfo[7].find("SRF=");
            size_t posSRFend = lineinfo[7].find(";", posSRFsta+1);
            if(posSRFend == std::string::npos)
            {
                posSRFend = lineinfo[7].size()+1;
            }
            string SRFCov    = lineinfo[7].substr(posSRFsta+4, posSRFend - posSRFsta - 4);
            double SRFcnt    = atof(SRFCov.c_str());
            //
            size_t posSRRsta = lineinfo[7].find("SRR=");
            size_t posSRRend = lineinfo[7].find(";", posSRRsta+1);
            if(posSRRend == std::string::npos)
            {
                posSRRend = lineinfo[7].size()+1;
            }
            string SRRCov    = lineinfo[7].substr(posSRRsta+4, posSRRend - posSRRsta - 4);  
            double SRRcnt    = atof(SRRCov.c_str());
            //
            size_t posSAFsta = lineinfo[7].find("SAF=");
            size_t posSAFend = lineinfo[7].find(";", posSAFsta+1);
            if(posSAFend == std::string::npos)
            {
                posSAFend = lineinfo[7].size()+1;
            }
            string SAFCov    = lineinfo[7].substr(posSAFsta+4, posSAFend - posSAFsta - 4);
            double SAFcnt    = atof(SAFCov.c_str());
            //
            size_t posSARsta = lineinfo[7].find("SAR=");
            size_t posSARend = lineinfo[7].find(";", posSARsta+1);
            if(posSARend == std::string::npos)
            {
                posSARend = lineinfo[7].size()+1;
            }
            string SARCov    = lineinfo[7].substr(posSARsta+4, posSARend - posSARsta - 4);
            double SARcnt    = atof(SARCov.c_str());
            maf              = (SAFcnt+SARcnt)/(SRFcnt+SRRcnt+SAFcnt+SARcnt);
            maf              = round(maf*1000000)/1000000;
        }
        else
        */ // the above is not used, as sometimes AB is 0 from longranger, but there are observed reads! why 20180629?
        if(lineinfo[7].find("RO=")!=std::string::npos && 
           lineinfo[7].find("AO=")!=std::string::npos)
        {
            size_t posROsta = lineinfo[7].find("RO=");
            size_t posROend = lineinfo[7].find(";", posROsta+1);
            if(posROend == std::string::npos)
            {
                posROend = lineinfo[7].size()+1;
            }
            string ROCov    = lineinfo[7].substr(posROsta+3, posROend - posROsta - 3);        
            size_t posAOsta = lineinfo[7].find("AO=");
            size_t posAOend = lineinfo[7].find(";", posAOsta+1);
            if(posAOend == std::string::npos)
            {
                posAOend = lineinfo[7].size()+1;
            }
            string AOCov    = lineinfo[7].substr(posAOsta+3, posAOend - posAOsta - 3);
            maf             = atof( AOCov.c_str() )/( atof( AOCov.c_str() ) + atof( ROCov.c_str() ) );
            maf             = round(maf*1000000)/1000000;
        }
        if(maf<minVarAF || maf>maxVarAF) // TODO: control with maximum as well! - 2018-03-12
        {
            /*
             cout << "   Check: skipped AF at " 
                  << lineinfo[0] << "\t" 
                  << lineinfo[1] << "\twith RO,AO:"
                  << ROCov       << ","             << AOCov << endl;     
            */
            continue; // caution: only vcf-variations with af>=0.2 used as markers.
        }
        // insert molecule info:
        // cluster barcodes in vcf into molecules
        vector<string>::iterator bcitr;
        vector<string>::iterator bcitr_end;
        // alternative from donor parent
        bcitr     = varBarcodes.begin();
        bcitr_end = varBarcodes.end();            
        while(bcitr != bcitr_end)
        {
            iMARKER tmpmkr;
            //tmpmkr.pos  = strtoul(lineinfo[1].c_str(), NULL, 0);
            tmpmkr.pos    = strtoul(dss.str().c_str(), NULL, 0);
            if(isGiven == true)
            {
                tmpmkr.from   = knownMarker[key];
            }
            else
            {
                tmpmkr.from   = "2";
            }
            tmpmkr.af     = maf;
            tmpmkr.cov    = posCov;
            tmpmkr.isGiven= isGiven;
            string bckey  = lineinfo[0] + "#" + *bcitr; // molecule id at chr: chr#bc
            (*moleculeSet).insert(std::pair<string, iMARKER>(bckey, tmpmkr));
            bcitr++;
        }
        // reference from another
        bcitr     = refBarcodes.begin();
        bcitr_end = refBarcodes.end();
        while(bcitr != bcitr_end)
        {
            iMARKER tmpmkr;
            //tmpmkr.pos  = strtoul(lineinfo[1].c_str(), NULL, 0);
            tmpmkr.pos    = strtoul(dss.str().c_str(), NULL, 0);
            if(isGiven == true && knownMarker[key].compare("1") == 0)
            {
                tmpmkr.from = "2";
            }
            else
            {
                tmpmkr.from = "1";
            }
            tmpmkr.af     = maf;
            tmpmkr.cov    = posCov;
            tmpmkr.isGiven= isGiven;                
            string bckey  = lineinfo[0] + "#" + *bcitr; // molecule id at chr: chr#bc
            (*moleculeSet).insert(std::pair<string, iMARKER>(bckey, tmpmkr));
            bcitr++;
        }
        // number of informative markers of different types intersecting given list
        if(isGiven == true)
        {
            if(itype.compare("del")==0) delnum ++;
            else
            if(itype.compare("ins")==0) insnum ++;
            else
            if(itype.compare("snp")==0) snpnum ++;
            else ;            
            mkrNumChr ++;
        }
    }
    cout << "   Info: number of collected markers    (excluding unknown) with cov>="  << minVarCov << "x: "  
         <<               mkrNumChr                           << endl;
    cout << "         "<< delnum        << " 1-bp deletions"  << endl;
    cout << "         "<< insnum        << " 1-bp insertions" << endl;
    cout << "         "<< snpnum        << "      snps"       << endl;
    cout << "   Info: number of checked   variations (including unknown) with cov>="  << minVarCov << "x: " 
         <<               mkrNumChr_all                       << endl;
    cout << "         "<< delnum_all    << " 1-bp deletions"  << endl;
    cout << "         "<< insnum_all    << " 1-bp insertions" << endl;
    cout << "         "<< snpnum_all    << "      snps"       << endl;    
    // close files
    fp.close();
    return true;
}
//
bool get_recombis_options(int                  argc, 
                         char*                 argv[],
                         string*               filevcf,
                         string*               fileMoleTable,
                         string*               fileReadNumDistr,
                         string*               fileMoleLenDistr,
                         string*               filemark1,
                         string*               filemark2,
                         int*                  min_marker,
                         int*                  min_span,
                         double*               min_score,
                         string*               outPrefix,
                         int*                  max_moleculeLen,
                         int*                  max_readnum,
                         int*                  max_CORDist,
                         int*                  min_CORDist,
                         int*                  max_COSize,
                         map<string, string>*  visitedVarChr,
                         map<string, string>*  visitedMolChr)
{
    int  ic;
    bool varf   = false;
    bool molf   = false;
    bool reaf   = false;
    bool lenf   = false;
    bool maf1   = false;
    bool maf2   = false;
    bool fsplit = false;
    int  n_chr  = 0;
    // log cmdline
    cout << "   CMDline: ";
    ic = 0;
    while(ic < argc)
    {
        cout << argv[ic] << " ";
	string optstr = (string)argv[ic];
	if(optstr.compare("--sma") == 0 || optstr.compare("--spa") == 0)
	{
	    int chr_n = get_chr_num_ii( (string)argv[ic+1] );
	    if( chr_n > n_chr )
	    {
		n_chr = chr_n;
	    }
	}
        ic ++;
    }
    cout << endl;
    if(n_chr == 0)
    {
	cout << "   Error: number of chrs detected in marker file is 0. " << endl;
	return false;
    }
    
    // get values
    ic = 2;    
    while (ic < argc) 
    {
        string optstr = (string)argv[ic];
        if(optstr.compare("--split") == 0)
        {
            ic ++;
            string tmpfileflag = (string)argv[ic];
            if(tmpflagstr.compare("tmptmp")==0)
            {
                vector<string> splitflaginfo = split_string(tmpfileflag, '/');
                tmpflagstr = splitflaginfo[splitflaginfo.size()-1];
            }
            string chrfilename(""); 
            // string tmpchr("12345");                // TODO: generalize to other number of chromosomes
	    map<int, string> chr_str;
	    for(int ichr = 0; ichr < n_chr; ichr ++)
	    {
	        std::stringstream ss;
		ss.str("");
	        ss << ichr + 1;
	        chr_str.insert(std::pair<int, string>(ichr, ss.str() ) );
	    }
            for(int ichr = 0; ichr < n_chr; ichr ++)
            {
                chrfilename = tmpfileflag + "_subfile_chr_" + chr_str[ichr] + "_longranger_variants.vcf.gz";
                igzstream iitfp;
                iitfp.open(chrfilename.c_str(), ios::in);
                if(iitfp.good())
                {
                    (*visitedVarChr).insert(std::pair<string, string>(chr_str[ichr], chrfilename));
                    iitfp.close();
                    cout << "   Info: with flag "
                         << tmpfileflag 
                         << ", found vcf file for chr " 
                       //<< tmpchr.substr(ichr, 1)
			 << chr_str[ichr]
                         << endl;
                }
                else
                {
                    cout << "   Warning: with flag " 
                         << tmpfileflag 
                         << ", cannot find vcf file for chr " 
                       //<< tmpchr.substr(ichr, 1)
			 << chr_str[ichr]
                         << endl;
                    //return false;
                }
                chrfilename.clear();
                chrfilename = tmpfileflag + "_subfile_chr_" + chr_str[ichr] + "_DrLink_molecules.txt.gz";
                iitfp.open(chrfilename.c_str(), ios::in);
                if(iitfp.good())
                {
                    (*visitedMolChr).insert(std::pair<string, string>(chr_str[ichr], chrfilename));
                    iitfp.close();
                    cout << "   Info: with flag "
                         << tmpfileflag 
                         << ", found molecule table file for chr " 
                       //<< tmpchr.substr(ichr, 1)
			 << chr_str[ichr]
                         << endl;                    
                }
                else
                {
                    cout << "   Warning: with flag " 
                         << tmpfileflag 
                         << " cannot find molecule table file for chr " 
                       //<< tmpchr.substr(ichr, 1)
			 << chr_str[ichr]
                         << endl;
                    //return false;
                }                
                chrfilename.clear();                
            }
            fsplit = true;
        }
        else
        if(optstr.compare("--var") == 0)
        {
            ic ++;
            *filevcf = (string)argv[ic];
            igzstream fp;
            fp.open((*filevcf).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: snp variant file provided: "    << *filevcf << endl;
                varf = true;
            }
            else
            {
                cout << "   Error: cannot open snp variant file " << *filevcf << endl;
                return false;
            }
        }
        else
        if(optstr.compare("--mol") == 0)
        {
            ic ++;
            *fileMoleTable = (string)argv[ic];
            igzstream fp;
            fp.open((*fileMoleTable).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: molecule meta info file provided: "    << *fileMoleTable << endl;
                molf = true;
            }
            else
            {
                cout << "   Error: cannot open molecule meta info file " << *fileMoleTable << endl;
                return false;
            }
        }
        else
        if(optstr.compare("--str") == 0)
        {
            ic ++;
            *fileReadNumDistr = (string)argv[ic];
            fstream fp;
            fp.open((*fileReadNumDistr).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: molecule read number file provided: "    << *fileReadNumDistr << endl;
                reaf = true;
            }
            else
            {
                cout << "   Error: cannot open molecule read number file " << *fileReadNumDistr << endl;
                return false;
            }
        }     
        else
        if(optstr.compare("--stl") == 0)
        {
            ic ++;
            *fileMoleLenDistr = (string)argv[ic];
            fstream fp;
            fp.open((*fileMoleLenDistr).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: molecule length file provided: "    << *fileMoleLenDistr << endl;
                lenf = true;
            }
            else
            {
                cout << "   Error: cannot open molecule length file " << *fileMoleLenDistr << endl;
                return false;
            }
        }
        else
        if(optstr.compare("--nur") == 0)
        {
            ic ++;
            *max_readnum = atoi(argv[ic]);
            if(*max_readnum>3)
            {
                cout << "   Info: maximum number of reads per molecule provided: "   << *max_readnum << endl;
            }
            else
            {
                cout << "   Error: maximum number of reads per molecule too small: " << *max_readnum << endl;
                return false;
            }
        }    
        else
        if(optstr.compare("--nul") == 0)
        {
            ic ++;
            *max_moleculeLen = atoi(argv[ic]);
            if(*max_moleculeLen>10000)
            {
                cout << "   Info: maximum length of molecules provided (bp): "   << *max_moleculeLen << endl;
            }
            else
            {
                cout << "   Error: maximum length of molecules too small (bp): " << *max_moleculeLen << endl;
                return false;
            }
        }    
        else
        if(optstr.compare("--acd") == 0)
        {
            ic ++;
            *max_CORDist = atoi(argv[ic]);
            if(*max_CORDist<100000)
            {
                cout << "   Info: maximum inter-read distance in CO provided (bp): "   << *max_CORDist << endl;
            }
            else
            {
                cout << "   Error: maximum inter-read distance in CO too large (bp): " << *max_CORDist << endl;
                return false;
            }
            // we'd be glad to see observed RD in CO to be smaller, that is, the interval is well covered by reads,
            // although these reads may not cover any markers.
        } 
        else
        if(optstr.compare("--icd") == 0)
        {
            ic ++;
            *min_CORDist = atoi(argv[ic]);
            if(*min_CORDist>=0)
            {
                cout << "   Info: minimum inter-read distance in CO provided (bp): "   << *min_CORDist << endl;
            }
            else
            {
                cout << "   Error: minimum inter-read distance in CO too small (bp): " << *min_CORDist << endl;
                return false;
            }
            // Given sufficient molecule base coverage, we'd be glad to see observed RD in CO to be smaller, 
            //    that is, the interval is well covered by reads, although these reads may not cover any markers. 
            // Given limited molecule base coverage (like 0.1-0.2x), if it's smaller than expected, this may indicate
            // false merged molecules.
        }  
        else
        if(optstr.compare("--aco") == 0)
        {
            ic ++;
            *max_COSize = atoi(argv[ic]);
            if(*max_COSize>0)
            {
                cout << "   Info: maximum size of CO provided (bp): "   << *max_COSize << endl;
            }
            else
            {
                cout << "   Error: maximum size of CO too small (bp): " << *max_COSize << endl;
                return false;
            }
        }                                 
        else
        if(optstr.compare("--spa") == 0)
        {
            ic ++;
            *filemark1 = (string)argv[ic];
            fstream fp;
            fp.open((*filemark1).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: parental marker file provided: "    << *filemark1 << endl;
                maf1 = true;
            }
            else
            {
                cout << "   Error: cannot open paternal marker file " << *filemark1 << endl;
                return false;
            }
        }
        else
        if(optstr.compare("--sma") == 0)
        {
            ic ++;
            *filemark2 = (string)argv[ic];
            fstream fp;
            fp.open((*filemark2).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: marental marker file provided: "    << *filemark2 << endl;
                maf2 = true;
            }
            else
            {
                cout << "   Error: cannot open maternal marker file " << *filemark2 << endl;
                return false;
            }
        } 
        else
        if(optstr.compare("--imn") == 0)
        {
            ic ++;
            *min_marker = atoi(argv[ic]);
            if(*min_marker>=2)
            {
                cout << "   Info: minimum number of allelic markers provided: "   << *min_marker << endl;
            }
            else
            {
                cout << "   Error: minimum number of allelic markers too small: " << *min_marker << endl;
                return false;
            }
        }  
        else
        if(optstr.compare("--imr") == 0)
        {
            ic ++;
            *min_span = atoi(argv[ic]);
            if(*min_span>=50)
            {
                cout << "   Info: minimum span of allelic markers provided: "   << *min_span << " bp" << endl;
            }
            else
            {
                cout << "   Error: minimum span of allelic markers too small: " << *min_span << " bp" << endl;
                return false;
            }
        }   
        else
        if(optstr.compare("--ims") == 0)
        {
            ic ++;
            *min_score = atof(argv[ic]);
            if(*min_score>=0.5)
            {
                cout << "   Info: minimum score of candidate co events provided: "   << *min_score << endl;
            }
            else
            {
                cout << "   Error: minimum score of candidate co events too small: " << *min_score << endl;
                return false;
            }
        } 
        else
        if(optstr.compare("--out") == 0)
        {
            ic ++;
            *outPrefix = string(argv[ic]);
            if((*outPrefix).size() > 0)
            {
                cout << "   Info: label of output files provided: " << *outPrefix << endl;
            }
            else
            {
                *outPrefix = "enjo_";
            }
        } 
        else
        if(optstr.compare("--ird") == 0)
        {
            // minimum number of read ids on each side of CO intercal, excluding the one on CO
            ic ++;
            minRead = atoi(argv[ic]);
            cout << "   Info: minimum number of uniq reads ids on each CO side (only when CO on the same fragment): "
                 << minRead << endl;
        }         
        else
        if(optstr.compare("--ivc") == 0)
        {
            // minimum coverage of variants when reading vcf
            ic ++;
            minVarCov = atoi(argv[ic]);
            cout << "   Info: minimum coverage of variant site used as markers provided: "   << minVarCov << endl;
        }  
        else
        if(optstr.compare("--avc") == 0)
        {
            // maximum coverage of variants when reading vcf
            ic ++;
            maxVarCov = atoi(argv[ic]);
            cout << "   Info: maximum coverage of variant site used as markers provided: "   << maxVarCov << endl;
        }          
        else
        if(optstr.compare("--iaf") == 0)
        {
            // minimum allele ferquency of varians when reading vcf
            ic ++;
            minVarAF = atof(argv[ic]);
            cout << "   Info: minimum allele frequency of variant site used as markers provided: "   << minVarAF << endl;
        } 
        else
        if(optstr.compare("--aaf") == 0)
        {
            // maximum allele ferquency of varians when reading vcf
            ic ++;
            maxVarAF = atof(argv[ic]);
            cout << "   Info: maximum allele frequency of variant site used as markers provided: "   << maxVarAF << endl;
        } 
        else
        if(optstr.compare("--rcm") == 0)
        {
            // to collect a marker, check if it is covered by an alignment
            rcm = true;
            cout << "   Info: checking if a marker is covered by an alignment is on." << endl;
        }                             
        else
        {
            cout << "   Warning: option " << argv[ic] << " was not recognized and ignored."  << endl;
        }
        // next option
        ic ++;
    }
    // check necessary files
    if(fsplit == false)
    {
        if(varf == false)
        {
            cout << "   Error: file missing, pls check --var. " << endl;
            return false;
        }
        if(molf == false)
        {
            cout << "   Error: file missing, pls check --mol. " << endl;
            return false;        
        } 
    }
    else
    {
        if((*visitedVarChr).size() != 0)
        {
            if(varf == true)
            {
                cout << "   Warning: both splitted and orginal vcf file provided with --var and --split. "
                     << endl
                     << "            splitted files will be used. "
                     << endl;
            }
        }
        else
        {
            if(varf == false)
            {
                cout << "   Warning: you need to set --var, as --split gave no vcf file. " << endl;
                return false;
            }            
        }
        if((*visitedMolChr).size() != 0)       
        {
            if(molf == true)
            {
                cout << "   Warning: both splitted and orginal mole file provided with --mol and --split. " 
                     << endl
                     << "            splitted files will be used. "
                     << endl;
            } 
        }
        else
        {
            if(molf == false)
            {
                cout << "   Warning: you need to set --mol, as --split gave no mole table file. " << endl;
                return false;
            }             
        }        
    }   
    if(reaf == false)
    {
        cout << "   Error: file missing, pls check --str. " << endl;
        return false;          
    }
    if(lenf == false)
    {
        cout << "   Error: file missing, pls check --stl. " << endl;
        return false;          
    }    
    if(maf1 == false)
    {
        cout << "   Error: file missing, pls check --spa. " << endl;
        return false;          
    }    
    if(maf2 == false)
    {
        cout << "   Error: file missing, pls check --sma. " << endl;
        return false;          
    }    
    return true;
}
int get_chr_num_ii(string marker_file)
{
    // here we get the number of chrs in the given marker list
    ifstream ifp;
    if(!ifp.good() )
    {
        cout << "   warning: file not found. " << endl;
        return 0;
    }
    map<string, int> chrs;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        vector<string> lineinfo = split_string(line, '\t');
        string chrid = lineinfo[1];
        if(chrs.find(chrid) ==chrs.end() )
        {
           chrs.insert(std::pair<string, int>(chrid, 1) );
        }
    }
    return chrs.size();
}
