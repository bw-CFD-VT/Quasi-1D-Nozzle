// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "Geometry.hpp"

using namespace std;

//---------------------------------------- Geometry Indexing -------------------------------------------------------------//
void Geometry_Indexing(double imax, double& dx, vector<double>& x_interface, vector<double>& x_cell_center)
{
    double Nozzle_Length = 2;
    double imin = -1;           //Start of converging section

    x_interface.resize(imax+1,0);
    x_cell_center.resize(imax,0); 
    
    dx = Nozzle_Length/imax; //Step size determined from number of cells selected
    double half_dx = dx/2;   

    for (double i = 0; i<imax+1; i++)
    {
        x_interface[i] = imin + i * dx;
    }


    for (int i = 0; i<imax; i++)
    {
        x_cell_center[i] = x_interface[i] + half_dx;
    }

    return;
}
//-----------------------------------------------------------------------------------------------------------------------//

//------------------------------------------ Nozzle Area at Cell Interface(s) -------------------------------------------// 
 void Area(vector<double> x, vector<double>& Area_x)
 {
     Area_x.resize(x.size(),0);

     for (int i = 0; i<x.size(); i++)
     {
         Area_x[i]=0.2+0.4*(1+sin(M_PI*(x[i]-0.5))); 
     }

    return;
 }
//-----------------------------------------------------------------------------------------------------------------------//