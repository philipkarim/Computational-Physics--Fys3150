#include "algoritmer.h"
#include"solsystem.h"
#include"celestialBody.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include<string>


using namespace  std;

int main()
{
  solsystem solsystem;

  int Numberofhs = 1000;
  double StartTime = 0.0;
  double FinalTime = 1.0;

  double vsc = (60 * 60 * 24 * 365) / (1.496e8); // Velocity scaling factor from km/s to au/year
  double x,y,vx,vy,m;
  x = 1.0; y = 0.0; vx = 0.0; vy = 2.0*pi; m = 6*pow(10,24);

  remove("Earth.txt");
  solsystem.lag_body(3e-6, 1, 0, 0, 2.0*pi, 2);  //Earth
  for(int i = 1;i < Numberofhs; i++ ){
    solsystem.kjoring_algoritme(solsystem.bodies[0],"Earth.txt", FinalTime, Numberofhs, 0);
  }

  return 0;
}
