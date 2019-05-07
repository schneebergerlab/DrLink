#include    <string>
#include <algorithm>
#include  <iostream>
using namespace std;

string check_allele_order(string leftAlleles, string rightAlleles)
{
    // order of alleles
    int left1  = std::count(leftAlleles.begin(), leftAlleles.end(), '1');
    int left2  = std::count(leftAlleles.begin(), leftAlleles.end(), '2');
    int right1 = std::count(rightAlleles.begin(), rightAlleles.end(), '1');
    int right2 = std::count(rightAlleles.begin(), rightAlleles.end(), '2');
    string ordergnome("");
    if(left1 > left2) ordergnome += "1";
    else
    if(left1 < left2) ordergnome += "2";
    else
    {
        ordergnome += "u";
        cout << "   Warning: equal alleles " << endl;
    }
    if(right1 > right2) ordergnome += "1";
    else
    if(right1 < right2) ordergnome += "2";
    else
    {
        ordergnome += "u";
        cout << "   Warning: equal alleles " << endl;
    } 
    if(ordergnome.compare("12")!=0 && ordergnome.compare("21")!=0)
    {
        // should not happen though!
        cout << "   Warning: incorrect order of alleles in prediction " 
             << leftAlleles << "\t" << rightAlleles << endl;
    }
    //
    return ordergnome;
}
