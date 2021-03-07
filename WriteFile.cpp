#include "WriteFile.hpp"

using namespace std;


void norm_file(string filename, int counter, int imax, vector<vector<double> > v)
{   
    ofstream file(filename,std::ios_base::app);
    if (file.is_open())
    {
        file<<counter<<",";
        //Write Variable Outputs to Corresponding Column/Variable Label
        for (int i = 0; i<3; i++)
        {   
            file << fixed << setprecision(14) << v[counter][i]<<",";
        }
            file << "\n";
            file.close();
    
    }
}

void mach_file(string filename, int counter, int imax, vector<vector<double> > v)
{   
    ofstream file(filename,std::ios_base::app);
    if (file.is_open())
    {
        file<<counter<<",";
        //Write Variable Outputs to Corresponding Column/Variable Label
        for (int i = 0; i<imax; i++)
        {   
            file << fixed << setprecision(14) << v[counter][i]<<",";
        }
            file << "\n";
            file.close();
    
    }
}