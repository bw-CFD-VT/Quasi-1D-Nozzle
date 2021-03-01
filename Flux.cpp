#include "Flux.hpp"

using namespace std;

void Flux (int counter,int NI, vector<vector<vector<double>>> V_cell_center, vector<vector<vector<double>>> &V_interface)
{

    V_interface.resize(counter+1);
    V_interface[counter].resize(NI,{0});

     for (int i = 1;i<NI-1;i++)
     { 
         V_interface[counter][i].resize(3);
         for (int j = 0; j<3; j++)
         {
            V_interface[counter][i][j] = 0.5*(V_cell_center[counter][i][j]+V_cell_center[counter][i-1][j]);
         }
            //  cout<<i<<", "<<V_interface[counter][i][0]<<", "<<V_interface[counter][i][1]<<", "<<V_interface[counter][i][2]<<endl;
     }

    return;
}

