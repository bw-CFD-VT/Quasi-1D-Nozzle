#include "Source_Term.hpp"

using namespace std;

void Source_Term (int counter,double imax,double dx, vector<double> Area_interface,
                  vector<vector<vector<double> > > V_cell_center,vector<vector<vector<double> > > &SourceTerm)
{
    SourceTerm.resize(counter+1);
    SourceTerm[counter].resize(imax,vector<double>(3,0));

    double delta_A_dx=0;

    for (int i = 0; i<imax; i++)
    {
        delta_A_dx = (Area_interface[i+1]-Area_interface[i])/dx;
        SourceTerm[counter][i][1] = V_cell_center[counter][i][2]*delta_A_dx;

        // cout<<SourceTerm[counter][i][0]<<", "<<SourceTerm[counter][i][1]<<", "<<SourceTerm[counter][i][2]<<endl;
    }

    return;
}