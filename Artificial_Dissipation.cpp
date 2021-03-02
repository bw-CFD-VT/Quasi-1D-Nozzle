#include "Artificial_Dissipation.hpp"
#include "Constants.hpp"    

using namespace std;

void Artifical_Dissipation (double K_2, double K_4,int counter,int imax, int NI, vector<vector<double>> lambda_max, 
                            vector<vector<vector<double>>> V_cell_center, vector<vector<vector<double>>> U_cell_center,
                            vector<vector<double>> V_ghost_inflow, vector<vector<double>> V_ghost_outflow,
                            vector<vector<vector<double>>> &d)
{

    d.resize(counter+1);
    d[counter].resize(NI);

    double lambda_half;
    double eps_half_2;
    double eps_half_4;

    // // Interior Nodes
    // for (int i = 2;i<NI-2;i++)
    //  { 
    //        d[counter][i].resize(3,0);
    //        lambda_half = lambda_max[counter][i+1] - lambda_max[counter][i-1]; 

    //        for (int j = 0; j<3; j++)

    //        {
             
    //          d[counter][i][j] = 

    //        }


    //  }








    return;
}
