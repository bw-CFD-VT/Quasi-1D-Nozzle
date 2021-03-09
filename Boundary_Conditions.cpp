// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "Boundary_Conditions.hpp"
#include "Isentropic_Flow.hpp"
#include "Constants.hpp"    

using namespace std;


/*
    Boundary Conditions determined for inflow and outflow using a ghost cell
    Inflow: 
            (1) Extrapolate Mach # from the interior and fix stagnation Temp. and Press.,
                use isentropic flow relations to determine primative variables
    Outflow: 
             (1) Supersonic: Extrapolate all primitive variables from internal cell center
             (2) Subsonic: Fix back pressure and extrapolate pressure and density from internal cell center
*/

void Boundary_Conditions(int counter,int imax, int Case_Flag,int ghost_cell,vector<vector<vector<double> > > V_cell_center,
                         vector<vector<double> > M_cell_center,vector<vector<vector<double> > > &V_Boundary,
                         vector<vector<double> > &M_Boundary,vector<vector<vector<double> > > &V_ghost_inflow,
                         vector<vector<vector<double> > > &V_ghost_outflow)
 {
    double p, rho, u, T;
    
    M_Boundary.resize(counter+1);
    M_Boundary[counter].resize(2,0);
    
    V_Boundary.resize(counter+1);
    V_Boundary[counter].resize(2,vector<double>(3,0));
     

     //---------------------------------------- INFLOW BC ------------------------------------------------------------// 
     vector<double> M_ghost_inflow(ghost_cell,0); 

     M_ghost_inflow[0] = 2*M_cell_center[counter][0]-M_cell_center[counter][1]; 

     if (M_ghost_inflow[0]<0)
     {
        M_ghost_inflow[0] = 0.0005;
     }
   
     M_Boundary[counter][0] = 0.5*(M_ghost_inflow[0]+M_cell_center[counter][0]);

   //   if (M_Boundary[counter][0]<0)
   //   {
   //      M_Boundary[counter][0] = 0.001;
   //   }

     Isentropic_Flow(M_Boundary[counter][0],rho, u, p, T);
     V_Boundary[counter][0] = {rho, u, p};

     V_ghost_inflow.resize(counter+1);
     V_ghost_inflow[counter].resize(ghost_cell,vector<double>(3,0));

     for (int i = 0; i<ghost_cell; i++)
     {
        Isentropic_Flow(M_ghost_inflow[i],rho, u, p, T);
        V_ghost_inflow[counter][i] = {rho, u, p};
     }
   //   cout<<V_ghost_inflow[counter][0][0]<<", "<<V_ghost_inflow[counter][0][1]<<", "<<V_ghost_inflow[counter][0][2]<<endl;
   //   cout<<M_Boundary[counter][0]<<endl;
   //   cout<<V_Boundary[counter][0][0]<<", "<<V_Boundary[counter][0][1]<<", "<<V_Boundary[counter][0][2]<<endl;
    //---------------------------------------------------------------------------------------------------------------//      

    //---------------------------------------- OUTFLOW BC -----------------------------------------------------------// 
     V_ghost_outflow.resize(counter+1);
     V_ghost_outflow[counter].resize(ghost_cell,vector<double>(3,0));

    //---------------------------------------- Supersonic Outflow BC -----------------------------------------------// 
     if (Case_Flag == 1) //Supersonic Outflow
     {
        for (int i = 0; i<3; i++)
        { 
            V_ghost_outflow[counter][0][i]=2*V_cell_center[counter][imax-1][i]-V_cell_center[counter][imax-2][i];
            V_Boundary[counter][1][i] = 0.5*(V_ghost_outflow[counter][0][i]+V_cell_center[counter][imax-1][i]);
        }
      //   cout<<V_ghost_outflow[counter][0][0]<<", "<<V_ghost_outflow[counter][0][1]<<", "<<V_ghost_outflow[counter][0][2]<<endl;
      //   cout<<V_Boundary[counter][1][0]<<", "<<V_Boundary[counter][1][1]<<", "<<V_Boundary[counter][1][2]<<endl;
        
     }

    //------------------------------------------- Subsonic Outflow BC -----------------------------------------------// 
    if (Case_Flag == 2) //Subsonic Outflow
     {
         
        for (int i = 0; i<3; i++)
        { 
            if (i==2)
            {
                V_Boundary[counter][1][i] = p_back;
                V_ghost_outflow[counter][0][i]=(2*p_back)-V_cell_center[counter][imax-1][i];
            }

            else
            {
                V_ghost_outflow[counter][0][i]=(2*V_cell_center[counter][imax-1][i])-V_cell_center[counter][imax-2][i];

                V_Boundary[counter][1][i] = 0.5*(V_ghost_outflow[counter][0][i]+V_cell_center[counter][imax-1][i]);
            }
            
        }
     }



    //---------------------------------------------------------------------------------------------------------------// 

    return;
 }
