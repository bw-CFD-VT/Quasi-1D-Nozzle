// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

// Solver Notes
//-------------------------------------//
// Central Flux Quadrature
// JST Artificial Dissipation based on implementation discussed in 
// 1 Ghost Cell at Inflow and Outflow

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
#include <stdio.h>

using namespace std;

int main()
{
    clear_existing_file(); //Delete existing files that will be written to

    int Case_Flag=1;

    int imax = 20;            //# of Cells --> ****EVEN # TO GET INTERFACE @ THROAT****
    


    //---------------------------------------- Nozzle Geometry Evaluation -----------------------------------------------------//
    vector<double> x_interface, x_cell_center,Area_interface,Area_cell_center;
    double dx;

    Geometry_Indexing(imax,dx,x_interface,x_cell_center);
    Area(x_interface,Area_interface);
    Area(x_cell_center,Area_cell_center);
    //-------------------------------------------------------------------------------------------------------------------------//
    
    //---------------------------------------- Exact Solution (Isentropic) ----------------------------------------------------//
    vector<vector<double> > V_exact(imax,vector<double>(3,0));
    vector<vector<double> > U_exact(imax,vector<double>(3,0));
    vector<double> M_exact(imax,0),rho_exact(imax,0), u_vel_exact(imax,0), p_exact(imax,0), T_exact(imax,0);
    Exact_Solution(imax,x_cell_center,Area_cell_center,M_exact,rho_exact,u_vel_exact,p_exact,T_exact,V_exact,U_exact);
    //-------------------------------------------------------------------------------------------------------------------------//


    //----------------------------------------------------Initial Conditions --------------------------------------------------//
    double ghost_cell = 1; // Initialize # of ghost cells to be used (refers to # per boundary, i.e. 1 = 1 @ inflow and 1 @ outflow)
    vector<vector<vector<double> > > V_cell_center(2,vector<vector<double> >(imax,vector<double>(3,0)));     // Matrix of Primative Variable Vectors at cell center, V = [V1,V2,V3]
    vector<vector<vector<double> > > U_cell_center(2,vector<vector<double> >(imax,vector<double>(3,0)));  // Matrix of Conserved Variable Vectors, U = [U1,U2,U3], [row = time][column = i]
    vector<double> M_cell_center(imax,0);              // Vector of Mach number at cell center 
    vector<vector<double> > V_Boundary(2,vector<double>(3,0));        // Vector Primative Variable Vectors at interface, V = [V1,V2,V3]
    vector<double> M_Boundary(2*ghost_cell,0);                 // Matrix of Mach number at interface, M = [M^n,Mn+1,M]
    vector<double> V_ghost_inflow(3,0);    // Matrix of Primative Variable Vectors @ inflow ghost cell(s), V = [V1,V2,V3]
    vector<double> V_ghost_outflow(3,0);   // Matrix of Primative Variable Vectors @ outfow ghost cell(s), V = [V1,V2,V3]
    vector<double> U_ghost_inflow(3,0);    // Matrix of Conserved Variable Vectors @ outfow ghost cell(s), U = [U1,U2,U3]
    vector<double> U_ghost_outflow(3,0);   // Matrix of Conserved Variable Vectors @ outfow ghost cell(s), U = [U1,U2,U3]

    Initial_Conditions(imax,x_cell_center,V_cell_center,M_cell_center);
    for (int i = 0; i<imax; i++) primative_to_conserved(V_cell_center[0][i],U_cell_center[0][i]);

    Boundary_Conditions(imax,Case_Flag,ghost_cell,V_cell_center,M_cell_center,V_Boundary,M_Boundary,V_ghost_inflow,V_ghost_outflow);
    primative_to_conserved(V_ghost_inflow,U_ghost_inflow);
    primative_to_conserved(V_ghost_outflow,U_ghost_outflow);
    //-------------------------------------------------------------------------------------------------------------------------//



    //------------------------------------------------- MAIN LOOP -------------------------------------------------------------//
    vector<vector<double> > F(imax+1,vector<double>(3,0));          // Vector of Flux Variables, F = [F1,F2,F3]
    vector<vector<double> > d(imax+1,vector<double>(3,0));          // Vector of Artificial Dissipation term, d = [d1,d2,d3] 
    vector<double> dt(imax,0);                                      // Vector of local step time
    vector<double> a(imax,0);                                       // Vector of local sound speed 
    vector<double> lambda_max(imax,0);                              // Vector of local max eigenvalue |u| + a
    vector<vector<double> > Residual(imax,vector<double>(3,0));     // Steady-State residual vector
    vector<vector<double> > SourceTerm(imax,vector<double>(3,0));   // Vectorof Source Term, S = [S1,S2,S3] [
    // double K_2 = 1.5*0.25;                                               // Aritifical Dissipation constant (2nd order damping term)
    // double K_4 = 0.015625;                                          // Aritifical Dissipation constant (4th order damping term)
    double K_2 = 0.01;                                               // Aritifical Dissipation constant (2nd order damping term)
    double K_4 = 1.5*0.015625; 

    //------------- Iterative Convergence Variable(s) ---------//
    int counter = 0;
    vector<double> L2_norm_residual(3,0);
    vector<double> L2_n (3,0);
    double convergence_criteria = 1e-10;
    // double CFL = .01;//SHOCK
    double CFL = .125;

    //----------------------------------------- MAIN LOOP -------------------------------------------------------------//  
    do
    {  
    if(counter%1000 == 0) cout<<"Iteration: "<<counter<<endl;
    Time_Step(imax,CFL,dx,V_cell_center,lambda_max,a,dt); 
    Flux(imax,V_Boundary,U_cell_center,F);
    Artifical_Dissipation(K_2,K_4,imax,ghost_cell,lambda_max,V_cell_center,U_cell_center,V_ghost_inflow,U_ghost_inflow,V_ghost_outflow,U_ghost_outflow,d);
    Source_Term(imax,dx, Area_interface,V_cell_center,SourceTerm);
    
    double d_t = *std::min_element(dt.begin(),dt.end()); //If global time step -> replace dt below

    //------------------- Explicit Euler Find U at step n +1 ---------------------------//
    for (int i = 0; i<imax; i++)
    {
        for (int j = 0; j<3; j++)
        {
            Residual[i][j] = ((F[i+1][j]+d[i+1][j])*Area_interface[i+1])
                                     -((F[i][j]+d[i][j])*Area_interface[i])
                                      -SourceTerm[i][j]*dx;

            U_cell_center[1][i][j] = U_cell_center[0][i][j] 
                                            -(d_t/(Area_cell_center[i]*dx))*Residual[i][j];
        }
    }
    //----------------------------------------------------------------------------------//    
    

   //---------------------- Calculate L2 Residual Norm --------------------------------//
    L2_Norm(imax,Residual,L2_norm_residual);
 
    if (counter == 0) //Initial Residual Norm Used to Scale Remaining Calculated Norms
    {
        L2_n[0] = L2_norm_residual[0];   L2_n[1] = L2_norm_residual[1];   L2_n[2] = L2_norm_residual[2];
    }

    for (int i = 0; i<3; i++)
    {
        L2_norm_residual[i] = L2_norm_residual[i]/L2_n[i];
    }
    //---------------------------------------------------------------------------------//


    //---------------------- Write Output ---------------------------------------------//
    if (counter%100 == 0)
    {
        residual_norm_file(residualnormfile,counter,imax,L2_norm_residual);
        mach_file(machfile,counter,imax,M_cell_center);
        prim_variable_file(rhofile,0,counter,imax,V_cell_center);
        prim_variable_file(ufile,1,counter,imax,V_cell_center);
        prim_variable_file(pressfile,2,counter,imax,V_cell_center);
        prim_variable_file(tempfile,3,counter,imax,V_cell_center);
    }
    //---------------------------------------------------------------------------------//


    //----------- Update Conserved + Primative Variable to be at step n ---------------//
    counter++; 
    U_cell_center[0] = U_cell_center[1];
    for (int i = 0; i<imax; i++) conserved_to_primative(U_cell_center[0][i],V_cell_center[0][i]);
    for (int i = 0; i<imax; i++) primative_to_conserved(V_cell_center[0][i],U_cell_center[0][i]); //Due to limiting prim. variables
    //---------------------------------------------------------------------------------//

    //---------------------- Update Mach Number ---------------------------------------//
    for (int i = 0; i<imax; i++)
    {
        double a_mach = 0;
        Sound_Speed(V_cell_center[0][i][0],V_cell_center[0][i][2],a_mach);
        M_cell_center[i] = V_cell_center[0][i][1]/a_mach;
    }
    //---------------------------------------------------------------------------------//

    //---------------------- Update Boundary Conditions -------------------------------//
    Boundary_Conditions(imax,Case_Flag,ghost_cell,V_cell_center,M_cell_center,V_Boundary,M_Boundary,V_ghost_inflow,V_ghost_outflow);
    primative_to_conserved(V_ghost_inflow,U_ghost_inflow);
    primative_to_conserved(V_ghost_outflow,U_ghost_outflow); 
    //---------------------------------------------------------------------------------//


    } while (L2_norm_residual[0]>convergence_criteria && L2_norm_residual[1] > convergence_criteria && L2_norm_residual[2] > convergence_criteria);
    //-------------------------------------------------------------------------------------------------------------------------//
    
    cout<<"Solution Converged at Iteration: "<< counter <<", "<<M_cell_center[imax-1]<<endl;

    //---------------------- Write Converged Solution(s)-------------------------------//
    residual_norm_file(residualnormfile,counter,imax,L2_norm_residual);
    mach_file(machfile,counter,imax,M_cell_center);
    prim_variable_file(rhofile,0,counter,imax,V_cell_center);
    prim_variable_file(ufile,1,counter,imax,V_cell_center);
    prim_variable_file(pressfile,2,counter,imax,V_cell_center);
    prim_variable_file(tempfile,3,counter,imax,V_cell_center);
    //---------------------------------------------------------------------------------//
    
   
    //-------------- Calculate DE for prim. + cons. variables (Case (1))---------------//
    vector<vector<double> > DE_U(imax,vector<double>(3,0));
    vector<vector<double> > DE_V(imax,vector<double>(3,0));
    vector<double> L2_norm_error(3,0);

    if (Case_Flag == 1)
    {

        for (int i = 0; i<imax; i++)
        {
            for (int j = 0; j<3; j++)
            {
                DE_U[i][j] = U_cell_center[0][i][j] - U_exact[i][j];
                DE_V[i][j] = V_cell_center[0][i][j] - V_exact[i][j];
            }

        }
        L2_Norm(imax,DE_U,L2_norm_error);
        error_norm_file(errornormfile_U,imax,L2_norm_error);

        L2_Norm(imax,DE_V,L2_norm_error);
        error_norm_file(errornormfile_V,imax,L2_norm_error);
    }
    //---------------------------------------------------------------------------------//
};