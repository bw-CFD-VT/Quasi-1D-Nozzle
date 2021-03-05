#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>
using namespace std;

void resize(int counter, vector<vector<vector<double> > > &x)
{
  x.resize(counter+1);
  x[counter].resize(2);

  for (int i = 0; i<2; i++)
  {
    x[counter][i].resize(2,0);
    for (int j = 0; j<2; j++)
    {
      x[counter][i][j] = counter+1;
    }
  }


  return;
}

int main()
{
    int counter = 0;

    vector<vector<vector<double> > > x;


    for (int i = 0; i<2; i++)
    {
      resize(counter,x);

      cout<<x[i][0][i]<<", "<<x[i][1][i]<<endl;
      counter++;
    }

    for (int i = 0; i<2; i++)
    {
      cout<<x[i][0][i]<<", "<<x[i][1][i]<<endl;
    }




}