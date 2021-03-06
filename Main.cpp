// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

// Central Flux Quadrature
// JST Artificial Dissipation

#include "Constants.hpp"
#include "Geometry.hpp"
#include "Initial_Conditions.hpp"
#include "Boundary_Conditions.hpp"
#include "Isentropic_Flow.hpp"
#include "TimeStep.hpp"
#include "SoundSpeed.hpp"
#include "VariableSwap.hpp"
#include "Flux.hpp"
#include "Artificial_Dissipation.hpp"
#include "Exact_Isentropic.hpp"
#include "Source_Term.hpp"
#include "L2_Norm.hpp"
#include "WriteFile.hpp"

using namespace std;

string filename = "Results.txt";

int main()
{
    int Case_Flag=1;
    // cout<<"\n"<<"Choose Case"<<"\n"<<"(1): Supersonic (Isentropic)"<<"\n"
    //        <<"(2): Subsonic (Shock in Nozzle)"<<"\n"<<endl;
    // cin>>Case_Flag;


    int imax = 10;            //# of Cells --> ****EVEN # TO GET INTERFACE @ THROAT****
    int NI = imax +1;        //Max # of Interfaces
    double dx;


    //---------------------------------------- Nozzle Geometry Evaluation -----------------------------------------------------//
    vector<double> x_interface, x_cell_center,Area_interface,Area_cell_center;

    Geometry_Indexing(imax,NI,dx,x_interface,x_cell_center);
    Area(x_interface,Area_interface);
    Area(x_cell_center,Area_cell_center);
    //-------------------------------------------------------------------------------------------------------------------------//

    //---------------------------------------- Exact Solution (Isentropic) ----------------------------------------------------//
    vector<double> M_exact(imax,0),rho_exact(imax,0), u_exact(imax,0), p_exact(imax,0), T_exact(imax,0);
    
    Isentropic_Nozzle_Exact(imax,x_cell_center,Area_cell_center,M_exact,rho_exact,u_exact,p_exact,T_exact);
    
    for (int i = 0; i<imax; i++)
    {
        cout<<M_exact[i]<<endl;
    }
    //-------------------------------------------------------------------------------------------------------------------------//

    

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

    double ghost_cell = 2; // Initialize # of ghost cells to be used (refers to # per boundary, i.e. 1 = 1 @ inflow and 1 @ outflow)
    int counter = 0;
    
    Initial_Conditions(imax,NI,x_cell_center,V_cell_center,M_cell_center);
    primative_to_conserved(counter,V_cell_center,U_cell_center);
    
    Boundary_Conditions(counter,Case_Flag,ghost_cell,imax,NI,V_cell_center,M_cell_center,V_Boundary,M_Boundary,V_ghost_inflow,V_ghost_outflow);
    primative_to_conserved(counter,V_ghost_inflow,U_ghost_inflow);
    primative_to_conserved(counter,V_ghost_outflow,U_ghost_outflow);
    //-------------------------------------------------------------------------------------------------------------------------//



    //------------------------------------------------- MAIN LOOP -------------------------------------------------------------//
  
    vector<vector<vector<double> > > F; // Matrix of Flux Variable Vectors, F = [F1,F2,F3], [row = time][column = i+1/2]
    vector<vector<vector<double> > > d; //Matrix of artificial dissipation vectors, d = [d1,d2,d3], [row = time][column = i+1/2]
    vector<vector<double> > dt; //Matrix of local step time, [row = time][column = i]
    vector<vector<double> > a;  //Matrix of local sound speed, [row = time][column = i]
    vector<vector<double> > lambda_max; //Matrix of local max eigenvalue |u| + a, [row = time][column = i]
    vector<vector<vector<double> > > Residual; //Steady-State residual vector, [row = time][column = i]
    vector<vector<vector<double> > > SourceTerm; //Matrix of Source Term, S = [S1,S2,S3] [row = time][column = i]

    //------------- Iterative Convergence Variable(s) ---------//
    vector<vector<double> > L2;
    vector<vector<vector<double> > > DE;
    vector<double> L2_n (3,0);

    double CFL = 0.1;
    double K_2 = 0.25;
    double K_4 = 0.015625;
    
      
    do
    {
        
    //---------------------- MAIN ITERATION ---------------------------------------//  
    Time_Step(counter,imax,CFL,dx,V_cell_center,lambda_max,a,dt); 
    Flux(counter,NI,V_Boundary,U_cell_center,F);
    Artifical_Dissipation(K_2,K_4,counter,imax,NI,ghost_cell,lambda_max,V_cell_center,U_cell_center,V_ghost_inflow,U_ghost_inflow,V_ghost_outflow,U_ghost_outflow,d);
    Source_Term(counter,imax,dx, Area_interface,V_cell_center,SourceTerm);
    
    Residual.resize(counter+1);
    Residual[counter].resize(imax,vector<double>(3,0));
    U_cell_center.resize(counter+2);
    U_cell_center[counter+1].resize(imax,vector<double>(3,0));
    double d_t = *std::min_element(dt[counter].begin(),dt[counter].end()); //If global time step -> replace dt below

    for (int i = 0; i<imax; i++)
    {
        
        for (int j = 0; j<3; j++)
        {
            Residual[counter][i][j] = (F[counter][i+1][j]+d[counter][i+1][j])*Area_interface[i+1]-
                                      (F[counter][i][j]+d[counter][i][j])*Area_interface[i]-SourceTerm[counter][i][j]*dx;

            U_cell_center[counter+1][i][j] = U_cell_center[counter][i][j] - (dt[counter][i]/(Area_cell_center[i]*dx))*Residual[counter][i][j];
        }
        //  cout<<Residual[counter][i][0]<<", "<<Residual[counter][i][1]<<", "<<Residual[counter][i][2]<<endl;
        //  cout<<U_cell_center[counter+1][i][0]<<", "<<U_cell_center[counter+1][i][1]<<", "<<U_cell_center[counter+1][i][2]<<endl;
    }
    //----------------------------------------------------------------------------//    
    

    //  if (counter > 2000) CFL = 0.1;
     if (counter > 2500)  CFL = 0.2;
     if (counter > 3000) CFL = 0.5;

    // cout<<"Calculating Norms"<<endl;
    L2_Norm(counter,imax,Residual,L2);
    
    if (counter == 0)
    {
        L2_n[0] = L2[counter][0];
        L2_n[1] = L2[counter][1];
        L2_n[2] = L2[counter][2];
    }
    
    for (int i = 0; i<3; i++)
    {
        // cout<<L2_n[i]<<endl;
        L2[counter][i] = L2[counter][i]/L2_n[i];
        // cout<<L2[counter][i]<<endl;
    }
    write_file(filename,counter,imax,L2);

    counter++;
    cout<<"Counter: "<<counter<<endl;

    conserved_to_primative(counter,U_cell_center,V_cell_center);

    //---------------------- Update Mach Number ---------------------------------------//
    M_cell_center.resize(counter+1);
    M_cell_center[counter].resize(imax,0);
    
    for (int i = 0; i<imax; i++)
    {
        double a_mach = 0;
        Sound_Speed(V_cell_center[counter][i][0],V_cell_center[counter][i][2],a_mach);
        M_cell_center[counter][i] = V_cell_center[counter][i][1]/a_mach;
         cout<<M_cell_center[counter][i]<<endl;
        // cout<<V_cell_center[counter][i][0]<<"\t"<<V_cell_center[counter][i][1]<<"\t"<<V_cell_center[counter][i][2]<<endl;
    }
    //---------------------------------------------------------------------------------//


    //---------------------- Update Boundary Conditions -------------------------------//
    Boundary_Conditions(counter,Case_Flag,ghost_cell,imax,NI,V_cell_center,M_cell_center,V_Boundary,M_Boundary,V_ghost_inflow,V_ghost_outflow);
    primative_to_conserved(counter,V_ghost_inflow,U_ghost_inflow);
    primative_to_conserved(counter,V_ghost_outflow,U_ghost_outflow); 
    //---------------------------------------------------------------------------------//
     
    } while (counter<5000);
    //-------------------------------------------------------------------------------------------------------------------------//
    cout<<"Broke Loop"<<endl;

};