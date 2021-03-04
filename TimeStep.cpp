#include "TimeStep.hpp"
#include "SoundSpeed.hpp"
#include "Constants.hpp"

using namespace std;

void Time_Step (int counter, int imax, double CFL, double dx, vector<vector<vector<double> > > V_cell_center,
                vector<vector<double> > &lambda_max, vector<vector<double> > &a,vector<vector<double> > &dt)
{

    a.resize(counter+1);
    a[counter].resize(imax,0);

    lambda_max.resize(counter+1);
    lambda_max[counter].resize(imax,0);

    dt.resize(counter+1);
    dt[counter].resize(imax,0);

    for (int i = 0; i<imax; i++)
    {
        Sound_Speed(V_cell_center[counter][i][0],V_cell_center[counter][i][2],a[counter][i]);
        lambda_max[counter][i] = abs(V_cell_center[counter][i][1])+ a[counter][i];
        dt[counter][i] = CFL*dx/lambda_max[counter][i];
       
    }

    return;
}