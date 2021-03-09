// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "L2_Norm.hpp"

using namespace std;

void L2_Norm(double imax,vector<vector<double> >Residual, vector<double> &L2)
{
    for (int j = 0; j<3; j++)
    {
        double sum = 0;
        for (int i = 0; i<imax; i++)
        {
           sum += Residual[i][j]*Residual[i][j];
        }
        L2[j] = sqrt(sum/imax);
    }

    return;

}