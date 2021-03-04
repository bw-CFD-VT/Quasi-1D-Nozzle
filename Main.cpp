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
#include "Source_Term.hpp"
#include "L2_Norm.hpp"

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

    Isentropic_Nozzle_Exact(NI,x_interface,Area_interface);

    //----------------------------------------------------Initial Conditions --------------------------------------------------//
    vector<vector<vector<double> > > V_cell_center;     // Matrix of Primative Variable Vectors at cell center, V = [V1,V2,V3]
    vector<vector<vector<double> > > U_cell_center; // Matrix of Conserved Variable Vectors, U = [U1,U2,U3], [row = time][column = i]
    vector<vector<double> > M_cell_center;              // Matrix of Mach number at cell center, M = [M^n,Mn+1,M]
    vector<vector<vector<double> > > V_Boundary;        // Matrix of Primative Variable Vectors at interface, V = [V1,V2,V3]
    vector<vector<double> > M_Boundary;                 // Matrix of Mach number at interface, M = [M^n,Mn+1,M]
    vector<vector<vector<double> > > V_ghost_inflow;    // Matrix of Primative Variable Vectors @ inflow ghost cell(s), V = [V1,V2,V3]
    vector<vector<vector<double> > > V_ghost_outflow;   // Matrix of Primative Variable Vectors @ outfow ghost cell(s), V = [V1,V2,V3]
    vector<vector<vector<double> > > U_ghost_inflow;    // Matrix of Conserved Variable Vectors @ outfow ghost cell(s), U = [U1,U2,U3]
    vector<vector<vector<double> > > U_ghost_outflow;   // Matrix of Conserved Variable Vectors @ outfow ghost cell(s), U = [U1,U2,U3]

    double ghost_cell = 1; // Initialize # of ghost cells to be used (refers to # per boundary, i.e. 1 = 1 @ inflow and 1 @ outflow)
    int counter = 0;
    
    Initial_Conditions(imax,NI,x_cell_center,V_cell_center,M_cell_center);
    primative_to_conserved(counter,V_cell_center,U_cell_center);
    
    Boundary_Conditions(counter,ghost_cell,imax,NI,V_cell_center,M_cell_center,V_Boundary,M_Boundary,V_ghost_inflow,V_ghost_outflow);
    primative_to_conserved(counter,V_ghost_inflow,U_ghost_inflow);
    primative_to_conserved(counter,V_ghost_outflow,U_ghost_outflow);


    //-------------------------------------------------------------------------------------------------------------------------//


    //---------------------------------------- MAIN LOOP --------------------------------------------------//
    vector<vector<vector<double> > > F; // Matrix of Flux Variable Vectors, F = [F1,F2,F3], [row = time][column = i+1/2]
    vector<vector<vector<double> > > d; //Matrix of artificial dissipation vectors, d = [d1,d2,d3], [row = time][column = i+1/2]
    vector<vector<double> > dt; //Matrix of local step time, [row = time][column = i]
    vector<vector<double> > a;  //Matrix of local sound speed, [row = time][column = i]
    vector<vector<double> > lambda_max; //Matrix of local max eigenvalue |u| + a, [row = time][column = i]
    vector<vector<vector<double> > > Residual; //Steady-State residual vector, [row = time][column = i]
    vector<vector<vector<double> > > SourceTerm; //Matrix of Source Term, S = [S1,S2,S3] [row = time][column = i]
    vector<vector<double> > L2;


    double CFL = 0.1;
    double K_2 = 1/2;
    double K_4 = 1/32;
    
      
    do
    {
        
    //-------------------------- Main Iteration -------------------------------//    
    Time_Step(counter,imax,CFL,dx,V_cell_center,lambda_max,a,dt); 
    Flux(counter,NI,V_Boundary,U_cell_center,F);
    // Artifical_Dissipation(K_2,K_4,counter,imax,NI,ghost_cell,lambda_max,V_cell_center,U_cell_center,V_ghost_inflow,U_ghost_inflow,V_ghost_outflow,U_ghost_outflow,d);
    Source_Term(counter,imax,dx, Area_interface,V_cell_center,SourceTerm);
    
    Residual.resize(counter+1);
    Residual[counter].resize(imax);
    U_cell_center.resize(counter+2);
    U_cell_center[counter+1].resize(imax);
 
    
    for (int i = 0; i<imax; i++)
    {
        Residual[counter][i].resize(3,0);
        U_cell_center[counter+1][i].resize(3,0);

        for (int j = 0; j<3; j++)
        {

            Residual[counter][i][j] = F[counter][i+1][j]*Area_interface[i+1]-F[counter][i][j]*Area_interface[i]-SourceTerm[counter][i][j]*dx;
            U_cell_center[counter+1][i][j] = U_cell_center[counter][i][j] - (dt[counter][i]/(Area_cell_center[i]*dx))*Residual[counter][i][j];
        
        }

    }
    cout<<"Calculating Norms"<<endl;
    L2_Norm(counter,imax,Residual,L2);

    counter++;

    conserved_to_primative(counter,U_cell_center,V_cell_center);

    M_cell_center.resize(counter+1);
    M_cell_center[counter].resize(imax,0);
    a.resize(counter+1);
    a[counter].resize(imax,0);

    for (int i = 0; i<imax; i++)
    {
        Sound_Speed(V_cell_center[counter][i][0],V_cell_center[counter][i][2],a[counter][i]);
        M_cell_center[counter][i] = V_cell_center[counter][i][1]/a[counter][i];
    }
    
    Boundary_Conditions(counter,ghost_cell,imax,NI,V_cell_center,M_cell_center,V_Boundary,M_Boundary,V_ghost_inflow,V_ghost_outflow);
    primative_to_conserved(counter,V_ghost_inflow,U_ghost_inflow);
    primative_to_conserved(counter,V_ghost_outflow,U_ghost_outflow); 

     cout<<"Counter: "<<counter<<endl;
    } while (counter < 100);
    cout<<"Broke Loop"<<endl;

};