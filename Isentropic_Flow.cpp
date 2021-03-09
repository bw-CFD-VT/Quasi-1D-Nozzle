// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "Isentropic_Flow.hpp"
#include "Constants.hpp"    

using namespace std;

void Isentropic_Flow (double Mach, double &rho, double &u, double &p, double &T)
{
    double psi = 1+((gam-1)/2)*Mach*Mach;
    T = T0/psi;                              // Static Temp,             T, Kelvin
    p = p0/pow(psi,(gam/(gam-1)));           // Static Press,            p, pa
    rho = p/(R*T);                           // Density,               rho, kg/m^3
    u = abs(Mach*sqrt(gam*R*T));             // x-comp. velocity,        u, m/s

    return;

}