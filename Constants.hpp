// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP
#include <string>

using namespace std;

//------------- Flow Problem Constants ---------------//
const double R_u = 8314.0;   //  J/(kmol*K)
const double MW_air = 28.96; //  kg/kmol
const double R = R_u/MW_air; //  J/kg-K
const double gam = 1.4;      // 
const double T0 = 600;       // K
const double p0 = 300e3;     // Pa
const double p_back = 120e3; // Pa
//---------------------------------------------------//

//--------------- File Output Names -----------------//
const string exactfile = "Exact_Solution_Isentropic";
const string residualnormfile = "residual_norm";
const string errornormfile_U = "error_norm_U";
const string errornormfile_V = "error_norm_V";
const string machfile = "mach";
const string rhofile = "rho";
const string ufile = "u";
const string pressfile = "press";
const string tempfile = "temp";
//---------------------------------------------------//

#endif