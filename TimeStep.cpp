// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "TimeStep.hpp"
#include "SoundSpeed.hpp"
#include "Constants.hpp"

using namespace std;

void Time_Step (int imax, double CFL, double dx, vector<vector<vector<double> > > V_cell_center,
                vector<double> &lambda_max, vector<double> &a,vector<double> &dt)
{

 
    for (int i = 0; i<imax; i++)
    {
        Sound_Speed(V_cell_center[0][i][0],V_cell_center[0][i][2],a[i]);
        lambda_max[i] = abs(V_cell_center[0][i][1])+ a[i];
        dt[i] = CFL*dx/lambda_max[i];
        // cout<<dt[i]<<endl;
        // cout<<lambda_max[i]<<endl;
       
    }

    return;
}