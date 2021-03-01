#include "TimeStep.hpp"
#include "SoundSpeed.hpp"
#include "Constants.hpp"

using namespace std;

void Time_Step (int counter, double CFL, double dx, vector<vector<vector<double>>> V_cell_center, vector<vector<double>> &dt)
{

    double a = 0; //sound speed
    double lambda_max = 0; //max eigenvalue |u| + a

    dt.resize(counter+1);
    dt[counter].resize(V_cell_center[counter].size());

    for (int i = 0; i<V_cell_center[counter].size(); i++)
    {
        Sound_Speed(V_cell_center[counter][i][0],V_cell_center[counter][i][2],a);
        // cout<<V_cell_center[counter][i][0]<<", "<< V_cell_center[counter][i][2]<<", "<<a<<endl;
        lambda_max = abs(V_cell_center[counter][i][1])+ a;
        // cout<<lambda_max<<", "<<counter<<endl;
        dt[counter][i] = CFL*dx/lambda_max;
        // cout<<dt[counter][i]<<endl;

       
       
    }

    return;
}