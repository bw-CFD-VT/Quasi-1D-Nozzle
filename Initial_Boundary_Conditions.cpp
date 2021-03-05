#include "Initial_Boundary_Conditions.hpp"
#include "Constants.hpp"    

using namespace std;

double rho, u, p, T;


/*
    Initial Conditions determine dusing a fixed Stagnation Pressure and Temperature

    Rho, T, p: Isentropic Relations
    Mach #: Per Section 3, Slide 8 --> (2) Linearly varying M

*/

void Initial_Conditions (int imax,int NI, vector<double> x_cell_center, vector<vector<vector<double> > > &V_cell_center, 
                        vector<vector<double> > &M_cell_center)
{

    M_cell_center.resize(1);
    M_cell_center[0].resize(imax,0);

    V_cell_center.resize(1);
    V_cell_center[0].resize(imax);
    


     for (int i = 0;i<imax;i++)
     {  
        M_cell_center[0][i] = 0.9*x_cell_center[i] + 1;                  //Linearly varying Mach #
        Isentropic_Relationships(M_cell_center[0][i],rho, u, p, T);
        V_cell_center[0][i] = {rho, u, p};                   // primative variable vector at location x at time t = 0
     }


    return;
}

/*
    Boundary Conditions determined for inflow and outflow using a ghost cell

    Inflow: 
            (1) Extrapolate Mach # from the interior and fix stagnation Temp. and Press.,
                use isentropic flow relations to determine primative variables


    Outflow: 
             (1) Supersonic: 
             (2) Subsonic: 
    Mach #: Per Section 3, Slide 8 --> (2) Linearly varying M

*/

void Boundary_Conditions(int counter,int Case_Flag,int ghost_cell,int imax,int NI,vector<vector<vector<double> > > V_cell_center,
                         vector<vector<double> > M_cell_center,vector<vector<vector<double> > > &V_Boundary,
                         vector<vector<double> > &M_Boundary,vector<vector<vector<double> > > &V_ghost_inflow,
                         vector<vector<vector<double> > > &V_ghost_outflow)
 {
     
    M_Boundary.resize(counter+1);
    M_Boundary[counter].resize(2,0);
    
    V_Boundary.resize(counter+1);
    V_Boundary[counter].resize(2,vector<double>(3,0));
     

     //---------------------------------------- INFLOW BC ------------------------------------------------------------// 
     vector<double> M_ghost_inflow(ghost_cell,0); 

     M_ghost_inflow[1] = 2*M_cell_center[counter][0]-M_cell_center[counter][1]; //Second ghost cell on LHS
     M_ghost_inflow[0] = 2*M_ghost_inflow[1]-M_cell_center[counter][0];

     if (M_ghost_inflow[1]<0)
     {
        M_ghost_inflow[0] = 0.01;
        M_ghost_inflow[1] = M_ghost_inflow[0] + 0.05;
     }
     if (M_ghost_inflow[0]<0)
     {
         M_ghost_inflow[0] = 0.001;
     }

     M_Boundary[counter][0] = 0.5*(M_ghost_inflow[1]+M_cell_center[counter][0]);
     Isentropic_Relationships(M_Boundary[counter][0],rho, u, p, T);
     V_Boundary[counter][0] = {rho, u, p};
     

     V_ghost_inflow.resize(counter+1);
     V_ghost_inflow[counter].resize(ghost_cell,vector<double>(3,0));
     for (int i = 0; i<ghost_cell; i++)
     {
        Isentropic_Relationships(M_ghost_inflow[i],rho, u, p, T);
        V_ghost_inflow[counter][i] = {rho, u, p};
        // cout<<M_ghost_inflow[i]<<endl;
        // cout<<V_ghost_inflow[counter][i][0]<<", "<<V_ghost_inflow[counter][i][1]<<", "<<V_ghost_inflow[counter][i][2]<<endl;
     }

    //  cout<<M_Boundary[counter][0]<<endl;
    //  cout<<V_Boundary[counter][0][0]<<", "<<V_Boundary[counter][0][1]<<", "<<V_Boundary[counter][0][2]<<endl;
    //---------------------------------------------------------------------------------------------------------------//      

    //---------------------------------------- OUTFLOW BC -----------------------------------------------------------// 
   
     V_ghost_outflow.resize(counter+1);
     V_ghost_outflow[counter].resize(ghost_cell,vector<double>(3,0));

    //---------------------------------------- Supersonic Outflow BC -----------------------------------------------// 
     if (Case_Flag == 1) //Supersonic Outflow
     {
        for (int i = 0; i<3; i++)
        { 
            V_ghost_outflow[counter][0][i]=(2*V_cell_center[counter][imax-1][i])-V_cell_center[counter][imax-2][i];
            V_Boundary[counter][1][i] = 0.5*(V_ghost_outflow[counter][0][i]+V_cell_center[counter][imax-1][i]);
            V_ghost_outflow[counter][1][i]=(2*V_ghost_outflow[counter][0][i]-V_cell_center[counter][imax-1][i]);
        }
     }
    //   cout<<V_Boundary[counter][1][0]<<", "<<V_Boundary[counter][1][1]<<", "<<V_Boundary[counter][1][2]<<endl;

    //------------------------------------------- Subsonic Outflow BC -----------------------------------------------// 
    // if (Case_Flag == 2) //Subsonic Outflow
    //  {
    //      double p_back = 120e3;
         
    //     for (int i = 0; i<3; i++)
    //     { 
    //         if (i==2)
    //         {
    //             V_Boundary[counter][1][i] = p_back;
    //             V_ghost_outflow[counter][0][i]=(2*p_back)-V_cell_center[counter][imax-1][i];

    //             if(ghost_cell==2)
    //             {
    //                 V_ghost_outflow[counter][1][i]=(2*V_ghost_outflow[counter][0][i]-V_cell_center[counter][imax-1][i]);
    //             }
    //         }

    //         else
    //         {
    //             V_ghost_outflow[counter][0][i]=(2*V_cell_center[counter][imax-1][i])-V_cell_center[counter][imax-2][i];

    //             if(ghost_cell==2)
    //             {
    //                 V_ghost_outflow[counter][1][i]=(2*V_ghost_outflow[counter][0][i]-V_cell_center[counter][imax-1][i]);
    //             }

    //             V_Boundary[counter][1][i] = 0.5*(V_ghost_outflow[counter][0][i]+V_cell_center[counter][imax-1][i]);
    //         }
            
    //     }
    //  }



    //---------------------------------------------------------------------------------------------------------------// 

    return;
 }



void Isentropic_Relationships (double Mach, double &rho, double &u, double &p, double &T)
{
    double psi = 1+((gam-1)/2)*Mach*Mach;
    T = T0/psi;                              // Static Temp,             T, Kelvin
    p = p0/pow(psi,(gam/(gam-1)));           // Static Press,            p, pa
    rho = p/(R*T);                           // Density,               rho, kg/m^3
    u = abs(Mach*sqrt(gam*R*T));             // x-comp. velocity,        u, m/s

    return;

}
