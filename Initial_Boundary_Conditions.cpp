#include "Initial_Boundary_Conditions.hpp"
#include "Constants.hpp"    

using namespace std;

double rho, u, p, T;


/*
    Initial Conditions determine dusing a fixed Stagnation Pressure and Temperature

    Rho, T, p: Isentropic Relations
    Mach #: Per Section 3, Slide 8 --> (2) Linearly varying M

*/

void Initial_Conditions (int imax,int NI,vector<double> x_cell_center, vector<double> x_interface,
                         vector<vector<vector<double>>> &V_cell_center, vector<vector<double>> &M_cell_center,
                         vector<vector<vector<double>>> &V_interface)
{

    M_cell_center.resize(1);
    M_cell_center[0].resize(imax,{0});

    V_cell_center.resize(1);
    V_cell_center[0].resize(imax,{0});
    


     for (int i = 0;i<imax;i++)
     {  
        M_cell_center[0][i] = 0.9*x_cell_center[i] + 1;                  //Linearly varying Mach #
        Isentropic_Relationships(M_cell_center[0][i],rho, u, p, T);
        V_cell_center[0][i] = {rho, u, p};                   // primative variable vector at location x at time t = 0

        // cout<<M_cell_center[0][i]<<endl;
        // cout<<V_cell_center[0][i][0]<<", "<<V_cell_center[0][i][1]<<", "<<V_cell_center[0][i][2]<<endl;
     }



    V_interface.resize(1);
    V_interface[0].resize(NI,{0});

     for (int i = 1;i<NI-1;i++)
     { 
         V_interface[0][i].resize(3);

         for (int j = 0; j<3; j++)
         {
            V_interface[0][i][j] = 0.5*(V_cell_center[0][i][j]+V_cell_center[0][i-1][j]);
         }
         //  cout<<i<<", "<<V_interface[0][i][0]<<", "<<V_interface[0][i][1]<<", "<<V_interface[0][i][2]<<endl;
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

void Boundary_Conditions(int ghost_cell,int imax,int NI,vector<double> x_cell_center, 
                         vector<double> x_interface,vector<vector<vector<double>>> V_cell_center,
                         vector<vector<double>> M_cell_center,vector<vector<vector<double>>> &V_interface,
                         vector<vector<double>> &M_interface)
 {

    M_interface.resize(1);
    M_interface[0].resize(NI,{0});

     //--------INFLOW-------//

     double M_ghost_inflow = 2*M_cell_center[0][0]-M_cell_center[0][1];
    //  cout<<M_ghost_inflow<<endl;

     if (M_ghost_inflow<0)
     {
         M_ghost_inflow = 0.01;
     }

     M_interface[0][0] = 0.5*(M_ghost_inflow+M_cell_center[0][0]);
     Isentropic_Relationships(M_interface[0][0],rho, u, p, T);
     V_interface[0][0] = {rho, u, p};
     
    //  cout<<M_interface[0][0]<<endl;
    //  cout<<V_interface[0][0][0]<<", "<<V_interface[0][0][1]<<", "<<V_interface[0][0][2]<<endl;
          



     //--------OUTFLOW -> SUPERSONIC-------//
     
     vector<double> V_ghost_outflow;
     V_ghost_outflow.resize(V_interface[0][0].size(),0);
     for (int i = 0; i<V_cell_center[0][0].size(); i++)
     {
        V_ghost_outflow[i]=2*V_cell_center[0][imax-1][i]-V_cell_center[0][imax-2][i];
        V_interface[0][NI-1][i] = 0.5*(V_ghost_outflow[i]+V_cell_center[0][imax-1][i]);

     }

    //   cout<<V_interface[0][NI-1][0]<<", "<<V_interface[0][NI-1][1]<<", "<<V_interface[0][NI-1][2]<<endl;




      //--------OUTFLOW -> SUBSONIC---------//


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
