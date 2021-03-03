#include "Artificial_Dissipation.hpp"
#include "VariableSwap.hpp"
#include "Constants.hpp"    

using namespace std;

void Artifical_Dissipation (double K_2, double K_4,int counter,int imax, int NI, int ghost_cell, vector<vector<double> > lambda_max, 
                            vector<vector<vector<double> > > V_cell_center, vector<vector<vector<double> > > U_cell_center,
                            vector<vector<vector<double> > > V_ghost_inflow, vector<vector<vector<double> > >U_ghost_inflow,
                            vector<vector<vector<double> > > V_ghost_outflow, vector<vector<vector<double> > >U_ghost_outflow,
                            vector<vector<vector<double> > > &d)
{


    //------------------- Setup Vector of Cell Pressure For Artificial Dissipation Calc --------------------//
    vector<double> Cell_Pressure(imax+(2*ghost_cell),0);


    for (int i = 0; i<ghost_cell; i++)
    {
        Cell_Pressure[i] = V_ghost_inflow[counter][i][2];
        Cell_Pressure[i+imax+ghost_cell] = V_ghost_outflow[counter][i][2];
        // cout<<V_ghost_inflow[counter][i][2]<<", "<<Cell_Pressure[i]<<"\n";
        // cout<<V_ghost_outflow[counter][i][2]<<", "<<Cell_Pressure[i+imax+2]<<"\n";
    }

    for (int i = ghost_cell; i<imax+ghost_cell; i++)
    {
        Cell_Pressure[i] = V_cell_center[counter][i-ghost_cell][2];
        // cout<<V_cell_center[counter][i-2][2]<<", "<<Cell_Pressure[i]<<"\n";
       
    }
    //------------------------------------------------------------------------------------------------------//

    //------------------- Setup Vector of Conserved Variables For Artificial Dissipation Calc --------------//
    vector<vector<double> > U_dissipation(imax+(2*ghost_cell));

    for (int i = 0; i<ghost_cell; i++)
    {
        for (int j = 0; j<3; j++)
        {
            U_dissipation[i].resize(3,0);
            U_dissipation[i+imax+ghost_cell].resize(3,0);
            U_dissipation[i][j] = U_ghost_inflow[counter][i][j];
            U_dissipation[i+imax+ghost_cell][j] = U_ghost_outflow[counter][i][j];
        }
    }
    
    for (int i = ghost_cell; i<imax+ghost_cell; i++)
    {
        for (int j = 0; j<3; j++)
        {
            U_dissipation[i].resize(3,0);
            U_dissipation[i][j] = U_cell_center[counter][i-ghost_cell][j];
        }
    }

    for (int i = 0; i<imax+(2*ghost_cell); i++)
    {
        cout<<U_dissipation[i][0]<<", "<<U_dissipation[i][1]<<", "<<U_dissipation[i][2]<<endl;
    }
    //------------------------------------------------------------------------------------------------------//


    d.resize(counter+1);
    d[counter].resize(NI);

    double lambda_half;
    double eps_half_2=0.001;
    double eps_half_4=0.002;


    for (int i = 1; i<NI-1; i++)
    {
        lambda_half = 0.5*(lambda_max[counter][i]+lambda_max[counter][i-1]);
        
        for (int j = 0; j<3; j++)
        {
            d[counter][i].resize(3,0);

            d[counter][i][j] = -(lambda_half*eps_half_2*(U_dissipation[i+1][j]-U_dissipation[i][j]) - 
                                (lambda_half*eps_half_4*(U_dissipation[i+2][j]-3*U_dissipation[i+1][j]+3*U_dissipation[i][j]-
                                 U_dissipation[i-1][j])));
        }

    }



    return;
}
