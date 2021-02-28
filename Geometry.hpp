#


#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>

using namespace std;

//---------------------------------------- Geometry Indexing -------------------------------------------------------------//
void Geometry_Indexing(double imax, double NI, double& dx, vector<double>& x_interface, vector<double>& x_cell_center);
//-----------------------------------------------------------------------------------------------------------------------//

//------------------------------------------ Nozzle Area at Cell Interface(s) -------------------------------------------// 
 void Area(vector<double> x, vector<double>& Area_x);
//-----------------------------------------------------------------------------------------------------------------------//
