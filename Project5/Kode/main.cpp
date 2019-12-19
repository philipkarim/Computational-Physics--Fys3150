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
int Numberofhs = 100;
double StartTime = 0.0;
double FinalTime = 1.0;
double vsc = (60 * 60 * 24 * 365) / (1.496e8); // Velocity scaling factor from km/s to au/year

int main()
{
  solsystem Solsystem;
  //
  // Solsystem.lag_body(1.6425e-7, 0.39, 0, 0, 48 * vsc, 2, 0); //Mercury
	// Solsystem.lag_body(2.4335e-6, 0.72, 0, 0, 35 * vsc,2,0); //Venus
	 Solsystem.lag_body(3e-6, 1, 0, 0, 2*pi, 2, 0);  //Earth
	// Solsystem.lag_body(3.195e-7, 1.524, 0, 0, 24.1*vsc,2,0); //Mars
	// Solsystem.lag_body(9.5e-4, 5.2, 0, 0, 13.07*vsc,2,0);  //Jupiter
	// Solsystem.lag_body(2.75e-4, 9.54, 0, 0, 9.6*vsc, 2,0); //Saturn
	// Solsystem.lag_body(4.4e-5, 19.19, 0, 0, 6.80*vsc, 2,0); //Uranus
	// Solsystem.lag_body(5.15e-5, 30.06, 0, 0, 5.43*vsc, 2,0); //Neptune
  remove("Energi.txt");
  remove("Mercury.txt");
  remove("Venus.txt");
  remove("Earth.txt");
  remove("Mars.txt");
  remove("Jupiter.txt");
  remove("Saturn.txt");
  remove("Uranus.txt");
  remove("Neptune.txt");


  for(int i = 1;i < Numberofhs; i++ ){
    // Solsystem.kjoring_algoritme(Solsystem.bodies[0], "Mercury.txt", FinalTime, Numberofhs, 0,0);
		// Solsystem.kjoring_algoritme(Solsystem.bodies[1], "Venus.txt", FinalTime, Numberofhs, 0,0);
		// Solsystem.kjoring_algoritme(Solsystem.bodies[2], "Earth.txt", FinalTime, Numberofhs, 0,0);
		 Solsystem.kjoring_algoritme(Solsystem.bodies[0], "Earth.txt", FinalTime, Numberofhs, 1,0);
     // Solsystem.kjoring_algoritme(Solsystem.bodies[3], "Earth.txt", FinalTime, Numberofhs, 0,0);

		// Solsystem.kjoring_algoritme(Solsystem.bodies[4], "Jupiter.txt", FinalTime, Numberofhs, 0,0);
		// Solsystem.kjoring_algoritme(Solsystem.bodies[5], "Saturn.txt", FinalTime, Numberofhs, 0,0);
		// Solsystem.kjoring_algoritme(Solsystem.bodies[6], "Uranus.txt", FinalTime, Numberofhs, 0,0);
		// Solsystem.kjoring_algoritme(Solsystem.bodies[7], "Neptune.txt", FinalTime, Numberofhs, 0,0);

  }
  //beregner og lagrer energi
  Solsystem.energi("Energi.txt",Numberofhs,FinalTime);

  // //beregner tiden
  // Solsystem.beregning_tid(Solsystem);
  //
  // //ser paa sol og merkur
  // solsystem sol_merkur;
  // sol_merkur.lag_body(1.6425e-7, 0.3075, 0, 0, 12.44, 2, 1);
  // sol_merkur.merkur_presesjon(Solsystem.bodies[0],100, 10000, "merkur_presesjon.txt");
  system("pause");
}
