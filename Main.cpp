// #include "FileOutput.hpp"
#include "Constants.hpp"
#include "Geometry.hpp"
#include "Initial_Boundary_Conditions.hpp"
#include "TimeStep.hpp"
#include "SoundSpeed.hpp"
#include "VariableSwap.hpp"
#include "Flux.hpp"

using namespace std;



int main()
{
   
    int imax = 6;            //# of Cells --> ****EVEN # TO GET INTERFACE @ THROAT****
    int NI = imax +1;        //Max # of Interfaces
    int t=0;                 // Initialize pseudo-time iteration counter/index
    double dx;


    //---------------------------------------- Nozzle Geometry Evaluation -----------------------------------------------------//
    vector<double> x_interface, x_cell_center,Area_interface,Area_cell_center;

    Geometry_Indexing(imax,NI,dx,x_interface,x_cell_center);
    Area(x_interface,Area_interface);
    Area(x_cell_center,Area_cell_center);
    //-------------------------------------------------------------------------------------------------------------------------//


    //----------------------------------------------------Initial Conditions --------------------------------------------------//
    vector<vector<vector<double>>> V_cell_center; // Matrix of Primative Variable Vectors at cell center, V = [V1,V2,V3]
    vector<vector<double>> M_cell_center;         // Matrix of Mach number at cell center, M = [M^n,Mn+1,M]
    vector<vector<vector<double>>> V_interface; // Matrix of Primative Variable Vectors at interface, V = [V1,V2,V3]
    vector<vector<double>> M_interface;         // Matrix of Mach number at interface, M = [M^n,Mn+1,M]
    double ghost_cell = 1; // Initialize # of ghost cells to be used (refers to # per boundary, i.e. 1 = 1 @ inflow and 1 @ outflow)
    int counter = 0;
    Initial_Conditions(imax,NI,x_cell_center,x_interface,V_cell_center,M_cell_center);
   
    Boundary_Conditions(counter,ghost_cell,imax,NI,x_cell_center,x_interface,V_cell_center,M_cell_center,V_interface,M_interface);
    
    //-------------------------------------------------------------------------------------------------------------------------//


    //---------------------------------------- MAIN LOOP --------------------------------------------------//
    
    vector<vector<vector<double>>> U; // Matrix of Conserved Variable Vectors, U = [U1,U2,U3]
    vector<vector<vector<double>>> F; // Matrix of Flux Variable Vectors, F = [F1,F2,F3]; 
    vector<vector<double>> dt;
    double CFL = 0.5;
    double a;
    do
    {
        
    Time_Step(counter,CFL,dx,V_cell_center,dt); 
    primative_to_conserved(counter,V_cell_center,U);
    Flux(counter, NI, V_cell_center,V_interface);

    // for (int i = 1;i<NI;i++)
    //  { 
    //        cout<<i<<", "<<V_interface[counter][i][0]<<", "<<V_interface[counter][i][1]<<", "<<V_interface[counter][i][2]<<endl;
    //  }


    counter++;
    for (int i = 0; i<imax; i++)
    {
        U[counter][0]
        // Sound_Speed(V_cell_center[counter][i][0],V_cell_center[counter][i][2],a);
        // M_cell_center[counter][i] = V_cell_center[counter][i][1]/a;
    }

    // Boundary_Conditions(counter,ghost_cell,imax,NI,x_cell_center,x_interface,V_cell_center,M_cell_center,V_interface,M_interface);
    
 

    // F.resize(counter+1);
    // F[counter].resize(NI);

    // for (int i = 0; i<NI; i++)
    // {
    //     F[counter][i]=0.5


    // }
    


     
     cout<<counter<<endl;
    } while (counter < 1);
    cout<<"Broke Loop"<<endl;

    //Flux
    //Resiudal Calc



};