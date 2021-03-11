// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "Artificial_Dissipation.hpp"
#include "VariableSwap.hpp"
#include "Constants.hpp"
#include <iomanip>    

using namespace std;

void Artifical_Dissipation (double K_2, double K_4,int imax, int ghost_cell, vector<double>lambda_max, 
                            vector<vector<vector<double> > > V_cell_center, vector<vector<vector<double> > > U_cell_center,
                            vector<double> V_ghost_inflow, vector<double> U_ghost_inflow,
                            vector<double> V_ghost_outflow, vector<double> U_ghost_outflow,
                            vector<vector<double> >&d)
{


    //------------------- Setup Vector of Cell Pressure For Artificial Dissipation Calc --------------------//
    vector<double> Cell_Pressure(imax+2*ghost_cell,0);

    Cell_Pressure[0] = V_ghost_inflow[2];
    Cell_Pressure[imax+ghost_cell] = V_ghost_outflow[2];


    for (int i = 0; i<imax; i++)
    {
        Cell_Pressure[i+ghost_cell] = V_cell_center[0][i][2];
    }
    //------------------------------------------------------------------------------------------------------//

    //------------------- Setup Vector of Conserved Variables For Artificial Dissipation Calc --------------//
    vector<vector<double> > U_dissipation(imax+2*ghost_cell,vector<double>(3,0));

    for (int i = 0; i<ghost_cell; i++)
    {
        for (int j = 0; j<3; j++)
        {
            U_dissipation[i][j] = U_ghost_inflow[j];
            U_dissipation[i+imax+ghost_cell][j] = U_ghost_outflow[j];
        }
    }
    
    for (int i = 0; i<imax; i++)
    {
        for (int j = 0; j<3; j++)
        {
            U_dissipation[i+ghost_cell][j] = U_cell_center[0][i][j];
        }
    }
    //------------------------------------------------------------------------------------------------------//

    //-------------------------------- Artificial Dissipation Term -----------------------------------------//
    double lambda_half;
    vector<double> v(imax+2*ghost_cell,0); // Pressure sensor vector

    double eps_half_2 = 0;
    double eps_half_4 = 0;

    for (int i = 1; i<imax+(2*ghost_cell)-1; i++)   
    {
      v[i] = abs(Cell_Pressure[i+1]-(2*Cell_Pressure[i])+Cell_Pressure[i-1])/
             abs(Cell_Pressure[i+1]+(2*Cell_Pressure[i])+Cell_Pressure[i-1]);
    }

    for (int i = 1; i<imax; i++)
    {   
        lambda_half = 0.5*(lambda_max[i]+lambda_max[i-1]);
        eps_half_2 = K_2*max_v(v[i-1],v[i],v[i+1],v[i+2]); 
        eps_half_4 = max(0.0,(K_4-eps_half_2));
      
        for (int j = 0; j<3; j++)
        {
            double D_1 = 0, D_3 = 0;
            D_1 = lambda_half*eps_half_2*(U_dissipation[i+1][j]-U_dissipation[i][j]);
            D_3= lambda_half*eps_half_4*(U_dissipation[i+2][j]-3*U_dissipation[i+1][j]+3*U_dissipation[i][j]-U_dissipation[i-1][j]);
            d[i][j] = D_3-D_1;
        }
    }

    for (int j = 0; j<3; j++)
    {
        d[0][j] = 2*d[1][j]-d[2][j];
        d[imax][j] = 2*d[imax-1][j]-d[imax-2][j];
    }

    return;
}

double max_v (double v_im1, double v_i, double v_ip1, double v_ip2)
{
    double max = v_im1;
    if(v_i>max) max=v_i;
    if(v_ip1>max) max=v_ip1;
    if(v_ip2>max) max=v_ip2;

    return (max);
}