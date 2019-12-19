#pragma once
#include "celestialBody.h"

#include <string>
#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <math.h>
#include<vector>


using namespace std;

class solsystem {
public:
  //definerer funksjoner 
  double const FourPi2 = 4 * pi*pi;
  vector<CelestialBody> bodies;
  solsystem();
  CelestialBody& lag_body(double masse, double xpos, double ypos, double xhas, double yhas, double beta,bool relativitisk_newton);
  void kalkulering_akselerasjon(vector<CelestialBody>& bodies);
  int numberOfBodies();
  void printBodies();
  void kjoring_algoritme(CelestialBody& body, string outfilename, double FinalTime, int Numberofhs,bool valg_av_algortime,bool relativitisk_newton);
  void merkur_presesjon (CelestialBody& body, double FinalTime, int Numberofhs,string outfilename);

};
