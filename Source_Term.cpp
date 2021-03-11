// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "Source_Term.hpp"

using namespace std;

void Source_Term (double imax,double dx, vector<double> Area_interface,
                  vector<vector<vector<double> > > V_cell_center,vector<vector<double> > &SourceTerm)
{

    double delta_A_dx=0;

    for (int i = 0; i<imax; i++)
    {
        delta_A_dx = (Area_interface[i+1]-Area_interface[i])/dx;
        SourceTerm[i][1] = V_cell_center[0][i][2]*delta_A_dx;
    }

    return;
}