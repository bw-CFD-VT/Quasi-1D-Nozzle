// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "Boundary_Conditions.hpp"
#include "Isentropic_Flow.hpp"
#include "Constants.hpp"    

using namespace std;


/*
    Boundary Conditions determined for inflow and outflow using 1 ghost cell
    Inflow: 
            (1) Extrapolate Mach # from the interior and fix stagnation Temp. and Press.,
                use isentropic flow relations to determine primitive variables
    Outflow: 
            (1) Supersonic: Extrapolate all primitive variables from internal cell center
            (2) Subsonic: Fix back pressure and extrapolate velocity and density from internal cell center
*/

void Boundary_Conditions(int imax, int Case_Flag,int ghost_cell,vector<vector<vector<double> > > V_cell_center,
                         vector<double> M_cell_center,vector<vector<double> > &V_Boundary,
                         vector<double> &M_Boundary,vector<double> &V_ghost_inflow, vector<double> &V_ghost_outflow)
 {
     double p, rho, u, T;
     
     //---------------------------------------- INFLOW BC ------------------------------------------------------------// 
     vector<double> M_ghost_inflow(ghost_cell,0); 

     M_ghost_inflow[0] = 2*M_cell_center[0]-M_cell_center[1]; 
     if (M_ghost_inflow[0]<0) M_ghost_inflow[0] = 0.0001;
     Isentropic_Flow(M_ghost_inflow[0],rho, u, p, T);
     V_ghost_inflow = {rho, u, p};

     M_Boundary[0] = 0.5*(M_ghost_inflow[0]+M_cell_center[0]);
     Isentropic_Flow(M_Boundary[0],rho, u, p, T);
     V_Boundary[0] = {rho, u, p};
    //---------------------------------------------------------------------------------------------------------------//    



    //---------------------------------------- OUTFLOW BC -----------------------------------------------------------// 

    //------------------------------- Case (1): Supersonic Outflow BC -----------------------------------// 
     if (Case_Flag == 1) //Supersonic Outflow
     {
        for (int i = 0; i<3; i++)
        { 
            V_ghost_outflow[i]=2*V_cell_center[0][imax-1][i]-V_cell_center[0][imax-2][i];
            V_Boundary[1][i] = 0.5*(V_ghost_outflow[i]+V_cell_center[0][imax-1][i]);
        }
     }
    //------------------------------ Case (2): Subsonic Outflow BC -------------------------------------// 
    if (Case_Flag == 2) //Subsonic Outflow
     {
         
        for (int i = 0; i<3; i++)
        { 
            if (i==2)
            {
                V_Boundary[1][i] = p_back;
                // V_ghost_outflow[i] = 2*V_Boundary[1][i]-V_cell_center[0][imax-1][i];
                V_ghost_outflow[i]=2*V_cell_center[0][imax-1][i]-V_cell_center[0][imax-2][i];
            }
            else
            {
                V_ghost_outflow[i]=(2*V_cell_center[0][imax-1][i])-V_cell_center[0][imax-2][i];
                V_Boundary[1][i] = 0.5*(V_ghost_outflow[i]+V_cell_center[0][imax-1][i]);
            }
        }
     }
    //---------------------------------------------------------------------------------------------------------------// 

    return;
 }
