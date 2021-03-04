#include "L2_Norm.hpp"

using namespace std;

void L2_Norm(int counter,double imax,vector<vector<vector<double> > > Residual, vector<vector<double> > &L2)
{

    L2.resize(counter+1);
    L2[counter].resize(3,0);


    for (int j = 0; j<3; j++)
    {
        double sum = 0;
        for (int i = 0; i<imax; i++)
        {
           sum += Residual[counter][i][j]*Residual[counter][i][j];
        }
        L2[counter][j] = sqrt(sum/imax);
    }

    return;

}