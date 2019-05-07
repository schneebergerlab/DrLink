#include <iostream>
#include  <fstream>
#include      <map>
#include   <string>
#include   <vector>

#include "split_string.h"

using namespace std;

bool read_chrsize(string chrSizeFile, map<string, string>* chrSizes)
{
    ifstream fp;
    fp.open(chrSizeFile.c_str(), ios::in);
    if(!fp.good())
    {
        cout << "   Error: cannot open chr size file " << chrSizeFile << endl;
        return false;
    }
    bool verbose = true;
    // Chr15	1091291
    cout << "\n   Info: reading sizes of chromosomes..." << endl;
    while(fp.good())
    {
        string line("");
        getline(fp, line);
        if(line.size()==0 || line[0]=='#') continue;
        
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 2) continue;
        
        map<string, string>::iterator itr;
        itr = (*chrSizes).find(lineinfo[0]);
        if(itr == (*chrSizes).end())
        {
            (*chrSizes).insert(std::pair<string, string>(lineinfo[0], lineinfo[1]));
            if(verbose)
            cout << "   Info: " << lineinfo[0] << " with size " << lineinfo[1]    << endl;
        }
        else
        {
            cout << "   Warning: repeated info " <<  line << endl;
        }        
    }
    fp.close();
    cout << "   Info: " << (*chrSizes).size() << " chrs recorded with size info." << endl << endl;
    return true;
}
