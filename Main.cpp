// #include "FileOutput.hpp"
#include "Geometry.hpp"
#include "Initial_Boundary_Conditions.hpp"
// #include "SoundSpeed.hpp"
#include "VariableSwap.hpp"
#include "Constants.hpp"

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


    //---------------------------------------- Boundary + Initial Conditions --------------------------------------------------//
    vector<vector<vector<double>>> V_cell_center; // Matrix of Primative Variable Vectors at cell center, V = [V1,V2,V3]
    vector<vector<double>> M_cell_center;         // Matrix of Mach number at cell center, M = [M^n,Mn+1,M]
    vector<vector<vector<double>>> V_interface; // Matrix of Primative Variable Vectors at interface, V = [V1,V2,V3]
    vector<vector<double>> M_interface;         // Matrix of Mach number at interface, M = [M^n,Mn+1,M]
    double ghost_cell = 1; // Initialize # of ghost cells to be used (refers to # per boundary, i.e. 1 = 1 @ inflow and 1 @ outflow)
    
    Initial_Conditions(imax,NI,x_cell_center,x_interface,V_cell_center,M_cell_center,V_interface);

    Boundary_Conditions(ghost_cell,imax,NI,x_cell_center,x_interface,V_cell_center,M_cell_center,V_interface,M_interface);
   
    for (int i = 0; i<NI; i++)
    {
        for (int j = 0; j<3; j++)
        {
            cout<<V_interface[0][i][j]<<", ";
        }
        cout<<"\n";
        
    }
    cout<<V_interface[0][1][0]<<endl;
   
    //-------------------------------------------------------------------------------------------------------------------------//

    
    
    // vector<vector<vector<double>>> U; // Matrix of Conserved Variable Vectors, U = [U1,U2,U3]
    // vector<vector<vector<double>>> F; // Matrix of Flux Variable Vectors, F = [F1,F2,F3]; 


    // TimeStep(dx, )
    //Time Step
    //Speed of Sound
    //Flux
    //Prim to Conserved.
    //Resiudal Calc



};