#include"celestialBody.h"
#include"solsystem.h"
#include "algoritmer.h"


#include <iostream>
#include <iomanip>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

CelestialBody::CelestialBody(double m, double x,double y,
  double vx,double vy,double b){
    masse = m;
    xpos = x;
    ypos = y;
    r = sqrt(xpos*xpos + ypos * ypos);
    xhas = vx;
    yhas = vy;
    beta = b;
    xaks = -(4 * pi*pi*xpos) / (pow(r, beta + 1));
		yaks = -(4 * pi*pi*ypos) / (pow(r, beta + 1));

  }
