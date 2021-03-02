// #include "FileOutput.hpp"
#include "Constants.hpp"
#include "Geometry.hpp"
#include "Initial_Boundary_Conditions.hpp"
#include "TimeStep.hpp"
#include "SoundSpeed.hpp"
#include "VariableSwap.hpp"
#include "Flux.hpp"
#include "Artificial_Dissipation.hpp"
#include "Exact_Isentropic.hpp"

using namespace std;



int main()
{
   
    int imax = 4;            //# of Cells --> ****EVEN # TO GET INTERFACE @ THROAT****
    int NI = imax +1;        //Max # of Interfaces
    int t=0;                 // Initialize pseudo-time iteration counter/index
    double dx;


    //---------------------------------------- Nozzle Geometry Evaluation -----------------------------------------------------//
    vector<double> x_interface, x_cell_center,Area_interface,Area_cell_center;

    Geometry_Indexing(imax,NI,dx,x_interface,x_cell_center);
    Area(x_interface,Area_interface);
    Area(x_cell_center,Area_cell_center);
    //-------------------------------------------------------------------------------------------------------------------------//

    Isentropic_Nozzle_Exact(NI,x_interface,Area_interface);

    //----------------------------------------------------Initial Conditions --------------------------------------------------//
    vector<vector<vector<double>>> V_cell_center; // Matrix of Primative Variable Vectors at cell center, V = [V1,V2,V3]
    vector<vector<double>> M_cell_center;         // Matrix of Mach number at cell center, M = [M^n,Mn+1,M]
    vector<vector<vector<double>>> V_Boundary; // Matrix of Primative Variable Vectors at interface, V = [V1,V2,V3]
    vector<vector<double>> M_Boundary;         // Matrix of Mach number at interface, M = [M^n,Mn+1,M]
    vector<vector<double>> V_ghost_inflow;
    vector<vector<double>> V_ghost_outflow;

    double ghost_cell = 1; // Initialize # of ghost cells to be used (refers to # per boundary, i.e. 1 = 1 @ inflow and 1 @ outflow)
    int counter = 0;
    
    Initial_Conditions(imax,NI,x_cell_center,V_cell_center,M_cell_center);
    Boundary_Conditions(counter,ghost_cell,imax,NI,V_cell_center,M_cell_center,V_Boundary,M_Boundary,V_ghost_inflow,V_ghost_outflow);
    
    //-------------------------------------------------------------------------------------------------------------------------//


    //---------------------------------------- MAIN LOOP --------------------------------------------------//
    
    vector<vector<vector<double>>> U_cell_center; // Matrix of Conserved Variable Vectors, U = [U1,U2,U3] [row = time][column = i]
    vector<vector<vector<double>>> F; // Matrix of Flux Variable Vectors, F = [F1,F2,F3] [row = time][column = i+1/2]
    vector<vector<vector<double>>> d; //Matrix of artificial dissipation vectors, d = [d1,d2,d3] [row = time][column = i+1/2]
    vector<vector<double>> dt; //Matrix of local step time [row = time][column = i]
    vector<vector<double>> a;  //Matrix of local sound speed [row = time][column = i]
    vector<vector<double>> lambda_max; //Matrix of local max eigenvalue |u| + a [row = time][column = i]

    double CFL = 0.5;
    double K_2 = 1/2;
    double K_4 = 1/32;
    
      
    do
    {
        
    //-------------------------- Main Iteration -------------------------------//    
    Time_Step(counter,imax,CFL,dx,V_cell_center,lambda_max,a,dt); 
    primative_to_conserved(counter,V_cell_center,U_cell_center);
    Flux(counter,NI,V_Boundary,U_cell_center,F);
    Artifical_Dissipation(K_2,K_4,counter,imax,NI,lambda_max,V_cell_center,U_cell_center,V_ghost_inflow,V_ghost_outflow,d);


    counter++;


     
     cout<<"Counter: "<<counter<<endl;
    } while (counter < 1);
    cout<<"Broke Loop"<<endl;

    //Flux
    //Resiudal Calc



};