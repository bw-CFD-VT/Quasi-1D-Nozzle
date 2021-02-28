#include "Geometry.hpp"

using namespace std;

//---------------------------------------- Geometry Indexing -------------------------------------------------------------//
void Geometry_Indexing(double imax, double NI, double& dx, vector<double>& x_interface, vector<double>& x_cell_center)
{
    double Nozzle_Length = 2;
    double imin = -1; //Min. x value of Nozzle, starting of converging section

    x_interface.resize(NI,0);
    x_cell_center.resize(imax,0); 

    dx = Nozzle_Length/imax; //Step size determined from number of "nodes" selected
    double half_dx = dx/2;


    for (double i = 0; i<NI; i++)
    {
        x_interface[i] = imin + i * dx;
        // cout<<x_interface[i]<<", ";
    }
        // cout<<"\n";

    for (double i = 0; i<imax; i++)
    {
        x_cell_center[i] = x_interface[i] + half_dx;
        // cout<<x_cell_center[i]<<", ";
    }
        // cout<<"\n";

    return;
}
//-----------------------------------------------------------------------------------------------------------------------//

//------------------------------------------ Nozzle Area at Cell Interface(s) -------------------------------------------// 
 void Area(vector<double> x, vector<double>& Area_x)
 {
     Area_x.resize(x.size());

     for (int i = 0; i<x.size(); i++)
     {
         Area_x[i]=0.2+0.4*(1+sin(M_PI*(x[i]-0.5))); 
        //  cout<<Area_x[i]<<", ";
     }
        // cout<<"\n";

    return;
 }
//-----------------------------------------------------------------------------------------------------------------------//
