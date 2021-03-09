// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "VariableSwap.hpp"
#include "Constants.hpp"

using namespace std;

//---------------------------------------- V -> U -----------------------------------------------------------------------//
void primative_to_conserved (int counter, vector<vector<vector<double> > > V, vector<vector<vector<double> > > &U)
{
    U.resize(counter+1);
    U[counter].resize(V[counter].size(),vector<double>(3,0));

        for (int j=0; j<V[counter].size(); j++)
        {
            U[counter][j][0] = V[counter][j][0];
            U[counter][j][1] = V[counter][j][0]*V[counter][j][1];
            U[counter][j][2] = (V[counter][j][2]/(gam-1))+(V[counter][j][0]*((V[counter][j][1]*V[counter][j][1])/2));
        }

    return;
}
//-----------------------------------------------------------------------------------------------------------------------//

//---------------------------------------- U -> V -----------------------------------------------------------------------//
void conserved_to_primative (int counter,vector<vector<vector<double> > > U, vector<vector<vector<double> > > &V)
{
    V.resize(counter+1);
    V[counter].resize(U[counter].size(),vector<double>(3,0));

        for (int j=0; j<U[counter].size(); j++)
        {
            V[counter][j][0] = U[counter][j][0];
            V[counter][j][1] = U[counter][j][1]/U[counter][j][0];
            V[counter][j][2] = (gam-1)*(U[counter][j][2] - U[counter][j][0]*(((U[counter][j][1]*U[counter][j][1])/(U[counter][j][0]*U[counter][j][0]))/2)); 
        }

    return;
}
//-----------------------------------------------------------------------------------------------------------------------//
