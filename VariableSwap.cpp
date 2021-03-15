// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "VariableSwap.hpp"
#include "Constants.hpp"

using namespace std;



//---------------------------------------- V -> U -----------------------------------------------------------------------//
void primitive_to_conserved (vector<double> V, vector<double> &U)
{
    U[0] = V[0];
    U[1] = V[0]*V[1];
    U[2] = (V[2]/(gam-1))+(V[0]*((V[1]*V[1])/2));

    return;
}
//-----------------------------------------------------------------------------------------------------------------------//

//---------------------------------------- U -> V -----------------------------------------------------------------------//
void conserved_to_primitive (vector<double> U, vector<double> &V)
{
    
    V[0] = U[0];
    V[1] = U[1]/U[0];
    V[2] = (gam-1)*(U[2] - U[0]*(((U[1]*U[1])/(U[0]*U[0]))/2)); 


    //------------ Limit Primative Variables ----------//
    vector<double> Primitive_Limits(3,0);
    
    Primitive_Limits[0] = 0.0001; Primitive_Limits[1] = 10.0; Primitive_Limits[2] = 500; //Nozzle Shock
    
    for (int j = 0; j<3; j++) V[j] = max(V[j],Primitive_Limits[j]);

    Primitive_Limits[0] = 5.5; Primitive_Limits[1] = 1000; Primitive_Limits[2] = 400e3; //Nozzle Shock

    for (int j = 0; j<3; j++) V[j] = min(V[j],Primitive_Limits[j]);
    //-------------------------------------------------//

    return;
}
//-----------------------------------------------------------------------------------------------------------------------//
