#include "SoundSpeed.hpp"
#include "Constants.hpp"

using namespace std;

void Sound_Speed (double rho, double p, double& a)
{
    a = sqrt(gam*p/rho);

    return;
}