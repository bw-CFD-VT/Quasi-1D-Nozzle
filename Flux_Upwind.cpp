// AOE 6145
// Homework 4: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "Constants.hpp"
#include "VariableSwap.hpp"
#include <algorithm> //std to include max functionality 

using namespace std;

void Flux_States(int imax, double eps, double kappa, vector<vector<double> > V_Boundary, vector<vector<vector<double> > > V_cell_center, 
                 vector<vector<double> > &U_interface_left, vector<vector<double> > &U_interface_right)
{
    vector<vector<double> > r_plus(imax+1,vector<double>(3,0));
    vector<vector<double> > r_minus(imax+1,vector<double>(3,0));
    vector<vector<double> > psi_plus(imax+1,vector<double>(3,0));
    vector<vector<double> > psi_minus(imax+1,vector<double>(3,0));
    vector<vector<double> > V_interface_left(imax+1,vector<double>(3,0));
    vector<vector<double> > V_interface_right(imax+1,vector<double>(3,0));
    vector<vector<double> > U_interface_left(imax+1,vector<double>(3,0));
    vector<vector<double> > U_interface_right(imax+1,vector<double>(3,0));


    double delta = 1e-6;

    for (int i = 1;i<imax+1;i++)
    {    
        for (int j = 0; j<3; j++)
        {   
            r_plus[i][j] = (V_cell_center[0][i+2][j]-V_cell_center[0][i+1][j])/(copysign(1.0,(V_cell_center[0][i+1][j]-V_cell_center[0][i][j]))*max(abs(V_cell_center[0][i+1][j]-V_cell_center[0][i][j]),delta));
            r_minus[i][j] = (V_cell_center[0][i][j]-V_cell_center[0][i-1][j])/(copysign(1.0,(V_cell_center[0][i+1][j]-V_cell_center[0][i][j]))*max(abs(V_cell_center[0][i+1][j]-V_cell_center[0][i][j]),delta));
        }
    }

    for (int i = 1;i<imax+1;i++)
    {    
        for (int j = 0; j<3; j++)
        {   
            psi_plus[i][j] = (r_plus[i][j]+abs(r_plus[i][j]))/(1+r_plus[i][j]);
            psi_minus[i][j] = (r_minus[i][j]+abs(r_minus[i][j]))/(1+r_minus[i][j]);
        }
    }

    for (int i = 1;i<imax+1;i++)
    {       
        for (int j = 0; j<3; j++)
        {
            V_interface_left[i][j] = V_cell_center[0][i][j]+(eps/4)*((1-kappa)*psi_plus[i][j]*(V_cell_center[0][i][j]-V_cell_center[0][i-1][j])
                                                                    +(1+kappa)*psi_minus[i][j]*(V_cell_center[0][i+1][j]-V_cell_center[0][i][j]));
        
            V_interface_right[i][j] = V_cell_center[0][i+1][j]-(eps/4)*((1-kappa)*psi_minus[i][j]*(V_cell_center[0][i+2][j]-V_cell_center[0][i+1][j])
                                                                    +(1+kappa)*psi_plus[i][j]*(V_cell_center[0][i+1][j]-V_cell_center[0][i][j]));
        }
    }

    for (int i = 0; i<imax; i++) primitive_to_conserved(V_interface_left[i],U_interface_left[i]);
    for (int i = 0; i<imax; i++) primitive_to_conserved(V_interface_right[i],U_interface_right[i]);

    return;
}

void Flux_VL()
{



    return;
}