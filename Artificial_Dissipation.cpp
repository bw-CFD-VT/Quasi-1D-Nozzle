#include "Artificial_Dissipation.hpp"
#include "VariableSwap.hpp"
#include "Constants.hpp"
#include <iomanip>    

using namespace std;

void Artifical_Dissipation (double K_2, double K_4,int counter,int imax, int NI, int ghost_cell, vector<vector<double> > lambda_max, 
                            vector<vector<vector<double> > > V_cell_center, vector<vector<vector<double> > > U_cell_center,
                            vector<vector<vector<double> > > V_ghost_inflow, vector<vector<vector<double> > >U_ghost_inflow,
                            vector<vector<vector<double> > > V_ghost_outflow, vector<vector<vector<double> > >U_ghost_outflow,
                            vector<vector<vector<double> > > &d)
{


    //------------------- Setup Vector of Cell Pressure For Artificial Dissipation Calc --------------------//
    vector<double> Cell_Pressure(imax+2*ghost_cell,0);

    for (int i = 0; i<ghost_cell; i++)
    {
        Cell_Pressure[i] = V_ghost_inflow[counter][i][2];
        Cell_Pressure[i+imax+ghost_cell] = V_ghost_outflow[counter][i][2];
    }

    for (int i = 0; i<imax; i++)
    {
        Cell_Pressure[i+ghost_cell] = V_cell_center[counter][i][2];
    }

    // for (int i = 0; i<imax+(2*ghost_cell); i++)
    // {
    //     cout<<Cell_Pressure[i]<<endl;
    // }
    //------------------------------------------------------------------------------------------------------//

    //------------------- Setup Vector of Conserved Variables For Artificial Dissipation Calc --------------//
    vector<vector<double> > U_dissipation(imax+2*ghost_cell,vector<double>(3,0));

    for (int i = 0; i<ghost_cell; i++)
    {
        for (int j = 0; j<3; j++)
        {
            U_dissipation[i][j] = U_ghost_inflow[counter][i][j];
            U_dissipation[i+imax+ghost_cell][j] = U_ghost_outflow[counter][i][j];
        }
    }
    
    for (int i = 0; i<imax; i++)
    {
        for (int j = 0; j<3; j++)
        {
            U_dissipation[i+ghost_cell][j] = U_cell_center[counter][i][j];
        }
    }

    // for (int i = 0; i<imax+(2*ghost_cell); i++)
    // {
    //     cout<<U_dissipation[i][0]<<"\t"<<U_dissipation[i][1]<<"\t"<<U_dissipation[i][2]<<"\t"<<endl;
    // }
    //------------------------------------------------------------------------------------------------------//


    d.resize(counter+1);
    d[counter].resize(NI,vector<double>(3,0));

    double lambda_half;
    vector<double> v(imax+ghost_cell,0); // Pressure sensor vector

    double eps_half_2;
    double eps_half_4;

    for (int i = 1; i<imax+(2*ghost_cell)-1; i++)   //Loop through from 2nd ghost cell @inflow to 2nd to last cell in domain
    {
      v[i-1] = abs(Cell_Pressure[i+1]-(2*Cell_Pressure[i])+Cell_Pressure[i-1])/
               abs(Cell_Pressure[i+1]+(2*Cell_Pressure[i])+Cell_Pressure[i-1]);
            // cout<<v[i-1]<<endl;
    }

    for (int i = 1; i<imax; i++)
    {   
        lambda_half = 0.5*(lambda_max[counter][i]+lambda_max[counter][i-1]);
        // cout<<lambda_half<<endl;
        eps_half_2 = K_2*max_v(v[i-1],v[i],v[i+1],v[i+2]);
        // cout<<eps_half_2<<endl; 
        eps_half_4 = max(0.0,(K_4-eps_half_2));
        // cout<<eps_half_4<<endl; 
        for (int j = 0; j<3; j++)
        {
            
            d[counter][i][j] = -((lambda_half*eps_half_2*(U_dissipation[i+2][j]-U_dissipation[i+1][j]))-
                                 (lambda_half*eps_half_4*(U_dissipation[i+3][j]-3*U_dissipation[i+2][j]+
                                  3*U_dissipation[i+1][j]-U_dissipation[i][j])));
        }
        //  cout<<d[counter][i][0]<<"\t"<<d[counter][i][1]<<"\t"<<d[counter][i][2]<<"\t"<<endl;
    }

    for (int j = 0; j<3; j++)
    {
        d[counter][0][j] = 2*d[counter][1][j]-d[counter][2][j];
        d[counter][imax][j] = 2*d[counter][imax-1][j]-d[counter][imax-2][j];
        //  cout<<d[counter][0][j]<<"\t"<<d[counter][imax][j]<<"\t"<<endl;
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