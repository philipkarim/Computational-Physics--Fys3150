#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include<string>
#include <chrono>

//Inkluderer headerfiles
#include "algoritmer.h"
#include"solsystem.h"
#include"celestialBody.h"


using namespace  std;

//Definerer steglengder, tider, og skaleringskonstant for hastighet
int Numberofhs = 10000;
double StartTime = 0.0;
double FinalTime = 200;
double konverteringsfaktor = (60 * 60 * 24 * 365) / (1.496e8); // denne faktoren konverterer fra km/s til AU/aar

int main()
{
  solsystem Solsystem;
  //her blir det lagt inn plantene som objekter
  Solsystem.lag_body(1.6425e-7, 0.39, 0, 0, 48 * konverteringsfaktor, 2, 0); //Merkur
	Solsystem.lag_body(2.4335e-6, 0.72, 0, 0, 35 * konverteringsfaktor,2,0); //venus
	Solsystem.lag_body(3e-6, 1, 0, 0, 2*pi, 2, 0);  //Jorda
	Solsystem.lag_body(3.195e-7, 1.524, 0, 0, 24.1*konverteringsfaktor,2,0); //Mars
	Solsystem.lag_body(9.5e-4, 5.2, 0, 0, 13.07*konverteringsfaktor,2,0);  //Jupiter
	Solsystem.lag_body(2.75e-4, 9.54, 0, 0, 9.6*konverteringsfaktor, 2,0); //Saturn
	Solsystem.lag_body(4.4e-5, 19.19, 0, 0, 6.80*konverteringsfaktor, 2,0); //Uranus
	Solsystem.lag_body(5.15e-5, 30.06, 0, 0, 5.43*konverteringsfaktor, 2,0); //Neptun

  //fjerner tidligere filer
  remove("Energi.txt");
  remove("merkur_fil.txt");
  remove("venus_fil.txt");
  remove("jorda_fil.txt");
  remove("mar_fil.txt");
  remove("jupiter_fil.txt");
  remove("saturn_fil.txt");
  remove("uranus_fil.txt");
  remove("neptun_fil.txt");

  //kjorer algoritmene med hensyn paa numberOfhs
  for(int i = 1;i < Numberofhs; i++ ){
    Solsystem.kjoring_algoritme(Solsystem.bodies[0], "merkur_fil.txt", FinalTime, Numberofhs, 0,0);
		Solsystem.kjoring_algoritme(Solsystem.bodies[1], "venus_fil.txt", FinalTime, Numberofhs, 0,0);
		Solsystem.kjoring_algoritme(Solsystem.bodies[2], "jorda_fil.txt", FinalTime, Numberofhs, 0,0);
    Solsystem.kjoring_algoritme(Solsystem.bodies[3], "mar_fil.txt", FinalTime, Numberofhs, 0,0);
		Solsystem.kjoring_algoritme(Solsystem.bodies[4], "jupiter_fil.txt", FinalTime, Numberofhs, 0,0);
		Solsystem.kjoring_algoritme(Solsystem.bodies[5], "saturn_fil.txt", FinalTime, Numberofhs, 0,0);
		Solsystem.kjoring_algoritme(Solsystem.bodies[6], "uranus_fil.txt", FinalTime, Numberofhs, 0,0);
		Solsystem.kjoring_algoritme(Solsystem.bodies[7], "neptun_fil.txt", FinalTime, Numberofhs, 0,0);

  }

  //ser paa sol og merkur
  solsystem sol_merkur;
  sol_merkur.lag_body(1.6425e-7, 0.3075, 0, 0, 12.44, 2, 1);
  sol_merkur.merkur_presesjon(Solsystem.bodies[0],100, 10000, "merkur_presesjon.txt");
  system("pause");
}
