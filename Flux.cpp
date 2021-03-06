#include "Flux.hpp"
#include "Constants.hpp"
#include "VariableSwap.hpp"

using namespace std;

void Flux (int counter,int NI, vector<vector<vector<double> > > V_Boundary,
           vector<vector<vector<double> > > U_cell_center, vector<vector<vector<double> > > &F)
{
    vector<vector<vector<double> > > U_interface;
    vector<vector<vector<double> > > U_boundary;

    U_interface.resize(counter+1);
    U_interface[counter].resize(NI,vector<double>(3,0));

    primative_to_conserved(counter,V_Boundary,U_boundary);
    U_interface[counter][0] = U_boundary[counter][0];
    U_interface[counter][NI-1] = U_boundary[counter][1];
    // cout<<U_interface[counter][0][0]<<"\t"<<U_interface[counter][0][1]<<"\t"<<U_interface[counter][0][2]<<endl;
    // cout<<U_interface[counter][NI-1][0]<<"\t"<<U_interface[counter][NI-1][1]<<"\t"<<U_interface[counter][NI-1][2]<<endl;

    for (int i = 1;i<NI-1;i++)
    {       
        for (int j = 0; j<3; j++)
        {
            U_interface[counter][i][j] = 0.5*(U_cell_center[counter][i][j]+U_cell_center[counter][i-1][j]);
        }
        // cout<<U_interface[counter][i][0]<<"\t"<<U_interface[counter][i][1]<<"\t"<<U_interface[counter][i][2]<<endl;
    }


    F.resize(counter+1);
    F[counter].resize(NI,vector<double>(3,0));

    for (int j=0; j<NI; j++)
    {
        F[counter][j][0] = U_interface[counter][j][1];
        F[counter][j][1] = ((3-gam)/2)*(U_interface[counter][j][1]*U_interface[counter][j][1]/U_interface[counter][j][0])+(gam-1)*U_interface[counter][j][2];
        F[counter][j][2] = U_interface[counter][j][2]*(U_interface[counter][j][1]/U_interface[counter][j][0])+(U_interface[counter][j][1]/U_interface[counter][j][0])*((gam-1)*U_interface[counter][j][2]-(((gam-1)/2)*(U_interface[counter][j][1]*U_interface[counter][j][1]/U_interface[counter][j][0])));
        // cout<<F[counter][j][0]<<"\t"<<F[counter][j][1]<<"\t"<<F[counter][j][2]<<endl;    
    }

    return;
}