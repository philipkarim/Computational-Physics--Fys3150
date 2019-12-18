#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include<string>
#include <chrono>


#include "algoritmer.h"
#include"solsystem.h"
#include"celestialBody.h"


using namespace std::chrono;
using namespace  std;

int main()
{
  solsystem Solsystem;

  int Numberofhs = 1000;
  double StartTime = 0.0;
  double FinalTime = 20;

  double vsc = (60 * 60 * 24 * 365) / (1.496e8); // Velocity scaling factor from km/s to au/year

  remove("Energi.txt");
  remove("Earth.txt");
  remove("Merkur.txt");

  Solsystem.lag_body(3e-6, 1, 0, 0, 2.0*pi, 2,0);  //Earth
  for(int i = 1;i < Numberofhs; i++ ){
    Solsystem.kjoring_algoritme(Solsystem.bodies[0],"Earth.txt", FinalTime, Numberofhs, 0);
  }
  //beregner og lagrer energi
  Solsystem.energi("Energi.txt",Numberofhs,FinalTime);

  //beregner tiden
  Solsystem.beregning_tid(Solsystem);

  ofstream myFile;
  myFile.open("Tid.txt");
	myFile << setprecision(8) << fixed << "N " << setw(20) << "Euler(s) " << setw(20) << "Verlet(s)" << endl;


		high_resolution_clock::time_point te1 = high_resolution_clock::now();
		Solsystem.kjoring_algoritme(Solsystem.bodies[0], "Earth.txt", FinalTime, 100000, 1);
		high_resolution_clock::time_point te2 = high_resolution_clock::now();
		duration<double> eulerTime = duration_cast<duration<double>>(te2 - te1);

		high_resolution_clock::time_point tv1 = high_resolution_clock::now();
		Solsystem.kjoring_algoritme(Solsystem.bodies[0], "Earth.txt", FinalTime, 100000, 0);
		high_resolution_clock::time_point tv2 = high_resolution_clock::now();
		duration<double> verletTime = duration_cast<duration<double>>(tv2 - tv1);

		myFile << setprecision(8) << fixed << 100000 << setw(20) << eulerTime.count() << setw(20) << verletTime.count() << endl;

  //ser paa sol og merkur
  solsystem sol_merkur;
  sol_merkur.lag_body(1.6425e-7, 0.3075, 0, 0, 12.44, 2, 1);
  sol_merkur.merkur_presesjon(Solsystem.bodies[0],100, 10000, "merkur_presesjon.txt");
  return 0;
}
