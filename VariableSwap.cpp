#include "VariableSwap.hpp"
#include "Constants.hpp"

using namespace std;

//---------------------------------------- V -> U -----------------------------------------------------------------------//
void primative_to_conserved (int counter, vector<vector<vector<double>>> V, vector<vector<vector<double>>> &U)
{
    U.resize(counter+1);
    for (int i=0; i<V.size(); i++)
    {
        U[i].resize(V[i].size(),{0});
        for (int j=0; j<V[i].size(); j++)
        {
            U[i][j].resize(V[i][j].size(),0);
            U[i][j][0] = V[i][j][0];
            U[i][j][1] = V[i][j][0]*V[i][j][1];
            U[i][j][2] = V[i][j][1]*(V[i][j][2]/(gam-1)+(V[i][j][0]*(V[i][j][1]*V[i][j][1])/2));
    
        }
    }
    return;
}
//-----------------------------------------------------------------------------------------------------------------------//

//---------------------------------------- U -> V -----------------------------------------------------------------------//
void conserved_to_primative (int counter,vector<vector<vector<double>>> U, vector<vector<vector<double>>> &V)
{
    V.resize(counter+1);
    
    for (int i=0; i<V.size(); i++)
    {
        V[i].resize(U[i].size(),{0});
        for (int j=0; j<U[i].size(); j++)
        {
            V[i][j].resize(U[i][j].size(),0);
            V[i][j][0] = U[i][j][0];
            V[i][j][1] = U[i][j][1]/U[i][j][0];
            V[i][j][2] = (gam-1)*((U[i][j][2]*U[i][j][0])/2 - (U[i][j][0]*U[i][j][1]*U[i][j][1]/(U[i][j][0]*U[i][j][0])/2));
            
        }
    }
    return;
}
//-----------------------------------------------------------------------------------------------------------------------//

// //---------------------------------------- U -> F (Central Quadrature)---------------------------------------------------//
// void conserved_to_primative (vector<vector<vector<double>>> U, vector<vector<vector<double>>> &F)
// {
    
//     for (int i=0; i<U.size(); i++)
//     {
//         V[i].resize(U[i].size(),{0});
//         for (int j=0; j<U[i].size(); j++)
//         {
//             V[i][j].resize(U[i][j].size(),0);
//             V[i][j][0] = U[i][j][0];
//             V[i][j][1] = U[i][j][1]/U[i][j][0];
//             V[i][j][2] = (gam-1)*((U[i][j][2]*U[i][j][0])/2 - (U[i][j][0]*U[i][j][1]*U[i][j][1]/(U[i][j][0]*U[i][j][0])/2));
            
//         }
//     }
//     return;
// }
// //-----------------------------------------------------------------------------------------------------------------------//