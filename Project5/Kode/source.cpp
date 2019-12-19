#include <cmath>
#include<string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace  std;
// output file as global variable
ofstream ofile;
// function declarations
void output( double, double, double, double, double);
void euler(double,double,double,double,double);

int main()
{
  int NumberofSteps = 1000;
  double FinalTime = 1.0;
  double Step = FinalTime/((double) NumberofSteps);
  double time = 0.0;
  double pi = acos(-1.0);
  double FourPi2 = 4*pi*pi;
  char *outfilename = "fil";

  ofile.open(outfilename);
  euler(time,FinalTime,Step, FourPi2,pi);
  ofile.close();
  return 0;
}
void euler(double time,double FinalTime,double Step, double FourPi2,double pi){
    //initial values
  double x =  1.0; double y =  0.0; double vx = 0.0; double vy = 2.0*pi;
  double r = sqrt(x*x+y*y);
  while (time <= FinalTime){
    x += Step*vx;
    y += Step*vy;
    vx -= Step*FourPi2*x/(r*r*r);
    vy -= Step*FourPi2*y/(r*r*r);
    r = sqrt(x*x+y*y);
    time += Step;
    output(time, x, y, vx, vy);   // write to file
  }
}

void velocity_verlet(){

}

void output(double time, double x, double y, double vx, double vy)
{
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << time;
  ofile << setw(15) << setprecision(8) << x;
  ofile << setw(15) << setprecision(8) << y;
  ofile << setw(15) << setprecision(8) << vx;
  ofile << setw(15) << setprecision(8) << vy << endl;
}
