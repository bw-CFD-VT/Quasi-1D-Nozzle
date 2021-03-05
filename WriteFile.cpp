#include "WriteFile.hpp"

using namespace std;


void write_file(string filename, int counter, int imax, vector<vector<vector<double> > > v)
{   
    ofstream file(filename,std::ios_base::app);
    if (file.is_open())
    {
        file<<counter<<"\n";
        //Write Variable Outputs to Corresponding Column/Variable Label
        for (int i = 0; i<imax; i++)
        {   
            for (int j = 0; j<3; j++)
            {
                file << fixed << setprecision(14) << v[counter][i][j]<<",";
            }
            file << "\n";
        }
            file.close();
    
    }
}