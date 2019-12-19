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
  double vx,double vy,double b,bool relativitisk_newton){
    //her registeres det variabler som skal brukes senere i solsystem
    masse = m;
    xpos = x;
    ypos = y;
    r = sqrt(xpos*xpos + ypos * ypos);
    xhas = vx;
    yhas = vy;
    beta = b;
    //beregnes akselerasjonen avhengig om det er relativitisk eller ikke
    if(!relativitisk_newton){
      xaks = -(4 * pi*pi*xpos) / (pow(r, beta + 1));
  		yaks = -(4 * pi*pi*ypos) / (pow(r, beta + 1));
    }else{
      xaks = -((4 * pi*pi*xpos) / (pow(r, beta + 1))) *(1+3*(pow(xpos*yhas,2) + pow(ypos*xhas,2)))/(pow(r,2)*pow(3e8*60*60*24*365/ 1.496e8, 2));
      yaks = -((4 * pi*pi*ypos) / (pow(r, beta + 1))) *(1+3*(pow(xpos*yhas,2) + pow(ypos*xhas,2)))/(pow(r,2)*pow(3e8*60*60*24*365/ 1.496e8, 2));
    }

  }
