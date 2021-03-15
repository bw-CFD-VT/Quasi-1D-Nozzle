// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "Flux.hpp"
#include "Constants.hpp"
#include "VariableSwap.hpp"

using namespace std;

void Flux (int imax, vector<vector<double> > V_Boundary, 
           vector<vector<vector<double> > > U_cell_center, vector<vector<double> > &F)
{
    vector<vector<double> > U_interface(imax+1,vector<double>(3,0));
    vector<vector<double> > U_Boundary(2,vector<double>(3,0));

    for(int i = 0; i<2; i++) primitive_to_conserved(V_Boundary[i],U_Boundary[i]);
    U_interface[0] = U_Boundary[0];
    U_interface[imax] = U_Boundary[1];


    for (int i = 1;i<imax;i++)
    {       
        for (int j = 0; j<3; j++)
        {
            U_interface[i][j] = 0.5*(U_cell_center[0][i][j]+U_cell_center[0][i-1][j]);
        }
    }

    for (int j=0; j<imax+1; j++)
    {
        F[j][0] = U_interface[j][1];
        F[j][1] = ((3-gam)/2)*(U_interface[j][1]*U_interface[j][1]/U_interface[j][0])+(gam-1)*U_interface[j][2];
        F[j][2] = U_interface[j][2]*(U_interface[j][1]/U_interface[j][0])+(U_interface[j][1]/U_interface[j][0])*((gam-1)*U_interface[j][2]-(((gam-1)/2)*(U_interface[j][1]*U_interface[j][1]/U_interface[j][0])));
    }

    return;
}