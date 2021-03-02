#include "Initial_Boundary_Conditions.hpp"
#include "Constants.hpp"    

using namespace std;

double rho, u, p, T;


/*
    Initial Conditions determine dusing a fixed Stagnation Pressure and Temperature

    Rho, T, p: Isentropic Relations
    Mach #: Per Section 3, Slide 8 --> (2) Linearly varying M

*/

void Initial_Conditions (int imax,int NI,vector<double> x_cell_center, vector<vector<vector<double>>> &V_cell_center, 
                        vector<vector<double>> &M_cell_center)
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

        // cout<<M_cell_center[0][i]<<endl;
        // cout<<V_cell_center[0][i][0]<<", "<<V_cell_center[0][i][1]<<", "<<V_cell_center[0][i][2]<<endl;
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

void Boundary_Conditions(int counter, int ghost_cell,int imax,int NI,vector<vector<vector<double>>> V_cell_center,
                         vector<vector<double>> M_cell_center,vector<vector<vector<double>>> &V_Boundary,
                         vector<vector<double>> &M_Boundary,vector<vector<double>> &V_ghost_inflow,
                         vector<vector<double>> &V_ghost_outflow)
 {

    M_Boundary.resize(counter+1);
    M_Boundary[counter].resize(2,0);

    V_Boundary.resize(counter+1);
    V_Boundary[counter].resize(2);
    V_Boundary[counter][0].resize(3,0); //Inflow Boundary 
    V_Boundary[counter][1].resize(3,0); //Outflow Boundary


     //--------INFLOW-------//

     double M_ghost_inflow = 2*M_cell_center[counter][0]-M_cell_center[counter][1];
    //  cout<<M_ghost_inflow<<endl;

     if (M_ghost_inflow<0)
     {
         M_ghost_inflow = 0.1;
     }

     M_Boundary[counter][0] = 0.5*(M_ghost_inflow+M_cell_center[counter][0]);
     Isentropic_Relationships(M_Boundary[0][0],rho, u, p, T);
     V_Boundary[counter][0] = {rho, u, p};
     

     V_ghost_inflow.resize(counter+1);
     V_ghost_inflow[counter].resize(3,0);
     Isentropic_Relationships(M_ghost_inflow,rho, u, p, T);
     V_ghost_inflow[counter] = {rho, u, p};

        // cout<<V_ghost_inflow[counter][0]<<", "<<V_ghost_inflow[counter][1]<<", "<<V_ghost_inflow[counter][2]<<endl;

    // cout<<M_Boundary[0][0]<<endl;
    //  cout<<V_Boundary[0][0][0]<<", "<<V_Boundary[0][0][1]<<", "<<V_Boundary[0][0][2]<<endl;
          

     //--------OUTFLOW -> SUPERSONIC-------//
     
   
     V_ghost_outflow.resize(counter+1);
     V_ghost_outflow[counter].resize(V_Boundary[counter][0].size());

     for (int i = 0; i<3; i++)
     {
        V_ghost_outflow[counter][i]=(2*V_cell_center[counter][imax-1][i])-V_cell_center[counter][imax-2][i];
        V_Boundary[counter][1][i] = 0.5*(V_ghost_outflow[counter][i]+V_cell_center[counter][imax-1][i]);
     }

    //    cout<<V_Boundary[0][1][0]<<", "<<V_Boundary[0][1][1]<<", "<<V_Boundary[0][1][2]<<endl;
    //  cout<<V_ghost_outflow[counter][0]<<", "<<V_ghost_outflow[counter][1]<<", "<<V_ghost_outflow[counter][2]<<endl;



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
