// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "Exact_Isentropic.hpp"
#include "Constants.hpp"
#include "iomanip"
using namespace std;

void Isentropic_Nozzle_Exact (double imax, vector<double> x, vector<double> Nozzle_Area, 
                              vector<double> &M_exact, vector<double> &rho_exact, vector<double> &u_exact,
                              vector<double> &p_exact, vector<double> &T_exact)
 {
     double es = 1e-10; 
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

         else               //Diverging Section -> Assuming Pb << Pe -> Supersonic
         {
             M_initial = 3;
             A_bar = Nozzle_Area[i]/A_star;
             Mach(M_initial,A_bar,es,M_final);
             M_exact[i] = M_final;
         }
        //  cout<<x[i]<<", "<<fixed<<setprecision(14)<<M[i]<<endl;

     }

     //-------------- Isentropic Nozzle Relationships --------------------//
     double psi=0;

     for (int i = 0;i<imax;i++)
     {
         psi = 1+((gam-1)/2)*M_exact[i]*M_exact[i];
         T_exact[i] = T0/psi;                         // Static Temp,    T, Kelvin
         p_exact[i] = p0/pow(psi,(gam/(gam-1)));      // Static Press,   p, pa
         rho_exact[i] = p_exact[i]/(R*T_exact[i]);                // Density,      rho, kg/m^3
         u_exact[i] = abs(M_exact[i]*sqrt(gam*R*T_exact[i]));     // x-velocity,     u, m/s
     }

     return;
 }



void Mach(double M_initial, double A_bar, double es, double &M)
{

    //Newton's Method Variable Definition
    double M_k = 0,         //Mach # at iteration k
           M_kp1 = 0,       //Mach # at iteration k+1
           phi_k = 0,       //Intermediate definition for area-mach # relation, reference lecture set 2 (exact solutions)
           F_k = 0,         //F(M) at iteration k
           F_prime_k = 0,   //partial F(m)/partial x at iteration k
           ea = 0;          //computed "error" at iteration k


    M_k = M_initial;    //Initialize Mach # at iteration k to be initial guess pending where in nozzle (supersonic vs. subsonic)
                        //Initial guess is passed from isentropic_nozzle function

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
        // cout<<F_k<<endl;
        counter++;                         //loop counter

            if (counter>10e3)  //Failsafe incase there was a reason it would get stuck iterating
            {
            cout<<"Stuck in Loop, Quitting";
            break;
            }

    }while(ea>es); //Once meets set convergence criteria quit loop

    M = M_kp1;

    return;  
}