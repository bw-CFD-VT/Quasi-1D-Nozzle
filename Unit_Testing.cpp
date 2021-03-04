/*

Unit Testing via Acutest: https://github.com/mity/acutest

Eliminates need for explicit testing framework i.e. cppunit

Strictly requires header function -> #include "acutest.hpp"

Tested Functions
    -Sound Speed
    -Time Step
    -Variable Swap -> Primative to Conserved | Conserved to Primative
    

*/


#include "acutest.hpp"
#include "Constants.hpp"
#include "Geometry.hpp"
#include "TimeStep.hpp"
#include "SoundSpeed.hpp"
#include "VariableSwap.hpp"
#include "Flux.hpp"
#include "Norm.hpp"
#include <iomanip>

 double error_tol = 1e-14;

void Test_Sound_Speed(void)
{
    double a;
    double rho = 1.225;
    double p = 101325.0;

    Sound_Speed(rho,p,a);
    
    double a_error = abs(340.29399054347107-a);
    TEST_ASSERT(a_error<error_tol);

}


void Test_Time_Step(void)
{
    int counter =0;
    double CFL = 0.5;
    double dx = 0.5;
    int imax = 1;

    vector<vector<vector<double> > > V;
    V.resize(counter+1);
    V[counter].resize(imax);

    for (double i = 0; i<imax; i++)
    {
        V[counter][i].resize(3);
        for (double j=0; j<3; j++)
        {
            V[counter][i][j] = j+1+i;
        }

    }
    vector<vector<double> > lambda_max;
    vector<vector<double> > a;
    vector<vector<double> > dt;
    Time_Step(counter,imax,CFL,dx,V,lambda_max,a,dt);

    
    //---------------- Verify lambda_max ---------------------//
    double lambda_max_error = abs(4.0493901531919200-lambda_max[counter][0]);
    TEST_ASSERT(lambda_max_error<error_tol);
    //--------------------------------------------------------//

    //---------------- Verify a (again) ----------------------//

    double a_error = abs(2.0493901531919200-a[counter][0]);
    TEST_ASSERT(a_error<error_tol);
    //--------------------------------------------------------//
   
    //---------------- Verify dt -----------------------------//
    double dt_error = abs(0.0617376914898996-a[counter][0]);
    TEST_ASSERT(a_error<error_tol);
    //--------------------------------------------------------//
    
}

void Test_Variable_Swap(void)
{
    int imax = 1;
    int counter = 0;
    vector<vector<vector<double> > > U;
    vector<vector<vector<double> > > V;
    V.resize(counter+1);
    V[counter].resize(imax);

    for (double i = 0; i<imax; i++)
    {
        V[counter][i].resize(3);
        for (double j=0; j<3; j++)
        {
            V[counter][i][j] = j+1.775;
            cout<<V[counter][i][j]<<endl;
        }

    }

    //---------------- Verify V -> U -------------------------//
    primative_to_conserved(counter,V,U);

    double U1_error= abs(1.7750000000000-U[counter][0][0]);
    double U2_error= abs(4.9256250000000-U[counter][0][1]);
    double U3_error= abs(16.2718046875000-U[counter][0][2]);

    TEST_ASSERT(U1_error<error_tol);
    TEST_ASSERT(U2_error<error_tol);
    TEST_ASSERT(U3_error<error_tol);
    //--------------------------------------------------------//

    //---------------- Verify U -> V -------------------------//
    conserved_to_primative(counter,U,V);

    double V1_error= abs(1.7750000000000-V[counter][0][0]);
    double V2_error= abs(2.7750000000000-V[counter][0][1]);
    double V3_error= abs(3.7750000000000-V[counter][0][2]);

    TEST_ASSERT(V1_error<error_tol);
    TEST_ASSERT(V2_error<error_tol);
    TEST_ASSERT(V3_error<error_tol);
    //--------------------------------------------------------//
   
}

TEST_LIST = { 
    {"Time Step",Test_Time_Step},
    {"Sound Speed",Test_Sound_Speed},
    {"Variable Swap",Test_Variable_Swap},
    {0}
    };