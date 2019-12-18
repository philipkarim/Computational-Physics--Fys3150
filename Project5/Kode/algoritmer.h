#pragma once
#include "celestialBody.h"
#include <fstream>
#include <iostream>
#include<string>

#include <iomanip>

void euler(CelestialBody& body, double FinalTime, int Numberofhs,string outfilename);
void velocityVerlet(CelestialBody& body, double FinalTime, int Numberofhs,string outfilename,bool relativitisk_newton);
void output( double, double, double, double, double,string,
  double,double, double ,double);
