// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "Exact_Isentropic.hpp"
#include "Constants.hpp"
#include "WriteFile.hpp"
#include "VariableSwap.hpp"
#include "Isentropic_Flow.hpp"
#include <iomanip>

using namespace std;


void Mach(double M_initial, double A_bar, double es, double &M)
{

    //Newton's Method Variable Definition
    double M_k = 0,         //Mach # at iteration k
           M_kp1 = 0,       //Mach # at iteration k+1
           phi_k = 0,       //Intermediate definition for area-mach # relation, reference lecture set 2 (exact solutions)
           F_k = 0,         //F(M) at iteration k
           F_prime_k = 0,   //partial F(m)/partial x at iteration k
           ea = 0;          //computed "error" at iteration k


    M_k = M_initial;    

    int counter = 0;         

    do
    {
        if (counter>0)
        {
          M_k = M_kp1;
        }
        
        phi_k = (2/(gam+1))*(1+((gam-1)/2)*M_k*M_k);
        F_k =  pow(phi_k,((gam+1)/(gam-1)))-(A_bar*A_bar*M_k*M_k);
        F_prime_k = 2*M_k*(pow(phi_k,(2/(gam-1)))-(A_bar*A_bar));
        
        M_kp1 = M_k - (F_k/F_prime_k);


        ea = abs((M_kp1-M_k)/M_kp1); //compute current iterations error
        // ea = F_k;                 //Wasn't working correctly -> no matter convergence tolerance M_exit ~ 3.22 not 3.16
        counter++;                         //loop counter

            if (counter>10e4)  //Failsafe incase there was a reason it would get stuck iterating
            {
            cout<<"Stuck in Loop, Quitting";
            break;
            }

    }while(ea>es); 

    M = M_kp1;

    return;  
}

void Isentropic_Nozzle_Exact (double imax, vector<double> x, vector<double> Nozzle_Area, 
                              vector<double> &M_exact, vector<double> &rho_exact, vector<double> &u_exact,
                              vector<double> &p_exact, vector<double> &T_exact)
 {
     double es = 1e-15; 
     double M_initial=0,A_bar=0,M_final=0; //Initialize variables required for Newton's method 
     
     double A_star = 0.2; //HARD CODED IN THROAT AREA BASED ON A(0.0)

     for (int i = 0; i<imax; i++)
     {
         if (x[i]<=0)     //Converging Section -> Subsonic
         {
             M_initial = 0.5;
             A_bar = Nozzle_Area[i]/A_star;
             Mach(M_initial,A_bar,es,M_final);
             M_exact[i] = M_final;
         }

         else            //Diverging Section -> Assuming Pb << Pe -> Supersonic
         {
             M_initial = 3;
             A_bar = Nozzle_Area[i]/A_star;
             Mach(M_initial,A_bar,es,M_final);
             M_exact[i] = M_final;
         }

     }

     //-------------- Isentropic Nozzle Relationships --------------------//
     for (int i = 0;i<imax;i++) Isentropic_Flow(M_exact[i],rho_exact[i],u_exact[i],p_exact[i],T_exact[i]);

     return;
 }


void Exact_Solution(int imax, vector<double> x_cell_center, vector<double> Area_cell_center,
                    vector<double> M_exact, vector<double> rho_exact, vector<double> u_vel_exact,
                    vector<double> p_exact, vector<double> T_exact, vector<vector<double> > V_exact, vector<vector<double> >&U_exact)
{
    
    Isentropic_Nozzle_Exact(imax,x_cell_center,Area_cell_center,M_exact,rho_exact,u_vel_exact,p_exact,T_exact);
    
    vector<vector<double> > exact_soln(7,vector<double>(imax,0));

    for (int i = 0; i<imax; i++)
    {
        exact_soln[0][i]=x_cell_center[i];
        exact_soln[1][i]=Area_cell_center[i];
        exact_soln[2][i]=M_exact[i];
        exact_soln[3][i]=rho_exact[i];
        exact_soln[4][i]=u_vel_exact[i];
        exact_soln[5][i]=p_exact[i];
        exact_soln[6][i]=T_exact[i];
        
        V_exact[i][0] = rho_exact[i];
        V_exact[i][1] = u_vel_exact[i];
        V_exact[i][2] = p_exact[i];
    } 

   for(int i = 0; i<imax; i++) primative_to_conserved(V_exact[i],U_exact[i]);

   exact_file(exactfile,exact_soln);
    
   return;
}


