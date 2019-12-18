#pragma once

#include <string>
#include <cmath>
using namespace std;
const double pi = 3.14159265359;

class CelestialBody {
public:
  //definerer variabler som skal brukes i cpp-filen

  double masse;
  double r;
  double xpos,ypos;
  double xhas,yhas;
  double xaks,yaks;
  double beta;

  CelestialBody(double m, double x,double y,double vx,double vy,double b,bool relativitisk_newton);

};
