// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "SoundSpeed.hpp"
#include "Constants.hpp"

using namespace std;

void Sound_Speed (double rho, double p, double& a)
{
    a = sqrt(gam*p/rho);

    return;
}