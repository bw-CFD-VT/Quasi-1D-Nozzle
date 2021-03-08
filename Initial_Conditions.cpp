#include "Initial_Conditions.hpp"
#include "Isentropic_Flow.hpp"
#include "Constants.hpp"    

using namespace std;

/*
    Initial Conditions determined using a fixed Stagnation Pressure and Temperature
    Rho, T, p: Isentropic Relations
    Mach #: Per Section 3, Slide 8 --> (2) Linearly varying M
*/

void Initial_Conditions (int imax, vector<double> x_cell_center, vector<vector<vector<double> > > &V_cell_center, 
                        vector<vector<double> > &M_cell_center)
{
    double rho, u, p, T;

    M_cell_center.resize(1);
    M_cell_center[0].resize(imax,0);

    V_cell_center.resize(1);
    V_cell_center[0].resize(imax);
    
     for (int i = 0;i<imax;i++)
     {  
        M_cell_center[0][i] = 0.9*x_cell_center[i] + 1;      // Linearly varying Mach #
        Isentropic_Flow(M_cell_center[0][i],rho, u, p, T);   // Isentropic Flow Assumed in Chamber/Plenum
        V_cell_center[0][i] = {rho, u, p};                   // Primative variable vector at cell center at time t = 0
     }


    return;
}