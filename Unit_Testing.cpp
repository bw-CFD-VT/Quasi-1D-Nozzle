// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

/*

Unit Testing via Acutest: https://github.com/mity/acutest

Eliminates need for explicit testing framework i.e. cppunit

Strictly requires header function -> #include "acutest.hpp"

"Unit Tested" Functions
    -Geometry -> Interface + Cell Center Location | Area at each corresponding x location
    -Sound Speed
    -Time Step
    -Variable Swap -> Primative to Conserved | Conserved to Primative
    -Norms
    -Flux

*/

#include "acutest.hpp"
#include "Constants.hpp"
#include "Geometry.hpp"
#include "TimeStep.hpp"
#include "SoundSpeed.hpp"
#include "VariableSwap.hpp"
#include "Flux.hpp"
#include "L2_Norm.hpp"
#include <iomanip>

double error_tol = 1e-14;

void Test_Geometry(void)
{
    int imax = 4;
    
    double dx;

    vector<double> x_interface(imax+1,0), x_cell_center(imax,0);
    vector<double> Area_interface(imax+1,0),Area_cell_center(imax,0);

    Geometry_Indexing(imax,dx,x_interface,x_cell_center);
    Area(x_interface,Area_interface);
    Area(x_cell_center,Area_cell_center);

    double dx_expected = 0.5000000;

    vector<double> x_cell_center_expected(imax,0),Area_cell_center_expected(imax,0),x_cell_center_error(imax,0),Area_cell_center_error(imax,0);
    x_cell_center_expected[0] = -0.75000;
    x_cell_center_expected[1] = -0.25000;
    x_cell_center_expected[2] = 0.25000;
    x_cell_center_expected[3] = 0.75000;
    Area_cell_center_expected[0] = 0.882842712474619;
    Area_cell_center_expected[1] = 0.317157287525381;
    Area_cell_center_expected[2] = 0.317157287525381;
    Area_cell_center_expected[3] = 0.882842712474619;

    vector<double> x_interface_expected(imax+1,0),Area_interface_expected(imax+1,0),x_interface_error(imax+1,0),Area_interface_error(imax+1,0);
    x_interface_expected[0] = -1.00000;
    x_interface_expected[1] = -0.50000;
    x_interface_expected[2] = 0.00000;
    x_interface_expected[3] = 0.50000;
    x_interface_expected[4] = 1.00000;
    Area_interface_expected[0] = 1.0000;
    Area_interface_expected[1] = 0.6000;
    Area_interface_expected[2] = 0.2000;
    Area_interface_expected[3] = 0.6000;
    Area_interface_expected[4] = 1.0000;
    
    for (int i = 0; i<imax; i++)
    {
        x_cell_center_error[i] = abs(x_cell_center[i]-x_cell_center_expected[i]);
        Area_cell_center_error[i] = abs(Area_cell_center[i]-Area_cell_center_expected[i]);
        TEST_ASSERT(x_cell_center_error[i]<error_tol);
        TEST_ASSERT(Area_cell_center_error[i]<error_tol);
    }
    for (int i = 0; i<imax+1; i++)
    {
        x_interface_error[i] = abs(x_interface[i]-x_interface_expected[i]);
        Area_interface_error[i] = abs(Area_interface[i]-Area_interface_expected[i]);
        TEST_ASSERT(x_interface_error[i]<error_tol);
        TEST_ASSERT(Area_interface_error[i]<error_tol);
    }

}

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
    double CFL = 0.5;
    double dx = 0.5;
    int imax = 1;

    vector<vector<vector<double> > > V(1,vector<vector<double> >(imax,vector<double>(3,0)));

    for (double i = 0; i<imax; i++)
    {
        for (double j=0; j<3; j++)
        {
            V[0][i][j] = j+1+i;
        }

    }
    vector<double> lambda_max(imax,0);
    vector<double> a(imax,0);
    vector<double> dt(imax,0);
    Time_Step(imax,CFL,dx,V,lambda_max,a,dt);

    
    //---------------- Verify lambda_max ---------------------//
    double lambda_max_error = abs(4.0493901531919200-lambda_max[0]);
    TEST_ASSERT(lambda_max_error<error_tol);
    //--------------------------------------------------------//

    //---------------- Verify a (again) ----------------------//

    double a_error = abs(2.0493901531919200-a[0]);
    TEST_ASSERT(a_error<error_tol);
    //--------------------------------------------------------//
   
    //---------------- Verify dt -----------------------------//
    double dt_error = abs(0.0617376914898996-dt[0]);
    TEST_ASSERT(dt_error<error_tol);
    //--------------------------------------------------------//
    
}

void Test_Variable_Swap_1(void) // Test for conserved and primitive variable vectors for cell center values (vector<vector<vector<>>>)
{
    int imax = 1;
    vector<vector<vector<double> > > U(1,vector<vector<double> >(imax,vector<double>(3,0)));
    vector<vector<vector<double> > > V(1,vector<vector<double> >(imax,vector<double>(3,0)));

    for (double i = 0; i<imax; i++)
    {
        for (double j=0; j<3; j++)
        {
            V[0][i][j] = j+1.775;
        }

    }

    //---------------- Verify V -> U -------------------------//
    primitive_to_conserved(V[0][0],U[0][0]);

    double U1_error= abs(1.7750000000000-U[0][0][0]);
    double U2_error= abs(4.9256250000000-U[0][0][1]);
    double U3_error= abs(16.2718046875000-U[0][0][2]);

    TEST_ASSERT(U1_error<error_tol);
    TEST_ASSERT(U2_error<error_tol);
    TEST_ASSERT(U3_error<error_tol);
    //--------------------------------------------------------//

    //---------------- Verify U -> V -------------------------//
    conserved_to_primitive(U[0][0],V[0][0]);

    double V1_error= abs(1.7750000000000-V[0][0][0]);
    double V2_error= abs(10.0000000000000-V[0][0][1]);
    double V3_error= abs(500.0000000000000-V[0][0][2]);

    TEST_ASSERT(V1_error<error_tol);
    TEST_ASSERT(V2_error<error_tol);
    TEST_ASSERT(V3_error<error_tol);
    //--------------------------------------------------------//
   
}

void Test_Variable_Swap_2(void) // Test for conserved and primitive variable vectors for boundary face values (<vector<vector<>>)
{
    int imax = 1;
    vector<vector<double> > U(imax,vector<double>(3,0));
    vector<vector<double> > V(imax,vector<double>(3,0));;

    for (double i = 0; i<imax; i++)
    {
        for (double j=0; j<3; j++)
        {
            V[0][j] = j+1.775;
        }

    }

    //---------------- Verify V -> U -------------------------//
    primitive_to_conserved(V[0],U[0]);

    double U1_error= abs(1.7750000000000-U[0][0]);
    double U2_error= abs(4.9256250000000-U[0][1]);
    double U3_error= abs(16.2718046875000-U[0][2]);

    TEST_ASSERT(U1_error<error_tol);
    TEST_ASSERT(U2_error<error_tol);
    TEST_ASSERT(U3_error<error_tol);
    //--------------------------------------------------------//

    //---------------- Verify U -> V -------------------------//
    conserved_to_primitive(U[0],V[0]);

    double V1_error= abs(1.7750000000000-V[0][0]);
    double V2_error= abs(10.0000000000000-V[0][1]);
    double V3_error= abs(500.0000000000000-V[0][2]);

    TEST_ASSERT(V1_error<error_tol);
    TEST_ASSERT(V2_error<error_tol);
    TEST_ASSERT(V3_error<error_tol);
    //--------------------------------------------------------//
   
}

void Test_Variable_Swap_3(void) // Test for conserved and primitive variable vectors for ghost cell values (<vector<>)
{
    int imax = 1;
    vector<double> U(3,0);
    vector<double> V(3,0);

    for (double i = 0; i<imax; i++)
    {
        for (double j=0; j<3; j++)
        {
            V[j] = j+1.775;
        }

    }

    //---------------- Verify V -> U -------------------------//
    primitive_to_conserved(V,U);

    double U1_error= abs(1.7750000000000-U[0]);
    double U2_error= abs(4.9256250000000-U[1]);
    double U3_error= abs(16.2718046875000-U[2]);

    TEST_ASSERT(U1_error<error_tol);
    TEST_ASSERT(U2_error<error_tol);
    TEST_ASSERT(U3_error<error_tol);
    //--------------------------------------------------------//

    //---------------- Verify U -> V -------------------------//
    conserved_to_primitive(U,V);

    double V1_error= abs(1.7750000000000-V[0]);
    double V2_error= abs(10.0000000000000-V[1]);
    double V3_error= abs(500.0000000000000-V[2]);

    TEST_ASSERT(V1_error<error_tol);
    TEST_ASSERT(V2_error<error_tol);
    TEST_ASSERT(V3_error<error_tol);
    //--------------------------------------------------------//
   
}

void Test_L2_Norm(void)
{
    int counter = 0;
    int imax = 5;
    vector<double> L2(imax,0);
    vector<vector<double> > Residual(imax,vector<double>(3,0));

    for (double i = 0; i<imax; i++)
    {
        for (double j=0; j<3; j++)
        {
            Residual[i][j] = j+1+i;
        }
    }


    L2_Norm(imax,Residual,L2);

    vector<double> L2_expected(3,0),L2_error(3,0);
    L2_expected[0] = 3.31662479035540;
    L2_expected[1] = 4.24264068711928;
    L2_expected[2] = 5.19615242270663;

    for (int i = 0;i<3;i++)
    {
        L2_error[i] = abs(L2_expected[i]-L2[i]);
        TEST_ASSERT(L2_error[i]<error_tol);
    }

}

void Test_Flux(void)
{
    int counter = 0;
    int imax = 3;

    vector<vector<double> > F(imax+1,vector<double>(3,0));
    vector<vector<double> > V(2,vector<double>(3,0));
    vector<vector<vector<double> > > U(1,vector<vector<double> >(imax,vector<double>(3,0)));

    for (double i = 0; i<2; i++)
    {
        for (double j=0; j<3; j++)
        {
            V[i][j] = j+1+i;
        }
    }

    for (double i = 0; i<imax; i++)
    {
        for (double j=0; j<3; j++)
        {
            U[0][i][j] = j+1+i;
        }
    }

    Flux(imax,V,U,F);

    vector<vector<double> > F_expected(imax+1,vector<double>(3,0));
    vector<vector<double> > F_error(imax+1,vector<double>(3,0));
    F_expected[0][0] = 2.000000000000000;
    F_expected[0][1] = 7.000000000000000;
    F_expected[0][2] = 25.000000000000000;
    F_expected[1][0] = 2.500000000000000;
    F_expected[1][1] = 4.733333333333330;
    F_expected[1][2] = 6.777777777777780;
    F_expected[2][0] = 3.500000000000000;
    F_expected[2][1] = 5.720000000000000;
    F_expected[2][2] = 7.448000000000000;
    F_expected[3][0] = 6.000000000000000;
    F_expected[3][1] = 22.000000000000000;
    F_expected[3][2] = 69.000000000000000;

    for (int i = 0; i<imax+1; i++)
    {
        for (int j = 0; j<3; j++)
        {
            F_error[i][j] = abs(F_expected[i][j] - F[i][j]);
            TEST_ASSERT(F_error[i][j]<error_tol);
        }
    }




}

TEST_LIST = { 
    {"Geometry",Test_Geometry},
    {"Time Step",Test_Time_Step},
    {"Sound Speed",Test_Sound_Speed},
    {"Variable Swap",Test_Variable_Swap_1},
    {"Variable Swap",Test_Variable_Swap_2},
    {"Variable Swap",Test_Variable_Swap_3},
    {"L2 Norm",Test_L2_Norm},
    {"Flux",Test_Flux},
    {0}
    };