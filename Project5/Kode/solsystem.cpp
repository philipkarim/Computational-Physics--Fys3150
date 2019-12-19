#include"solsystem.h"
#include"celestialBody.h"
#include "algoritmer.h"

#include <iostream>
#include <iomanip>
#include <math.h>
#include<string>

using namespace std;
using namespace std::chrono;

solsystem::solsystem(){

}

CelestialBody& solsystem::lag_body(double masse, double xpos, double ypos, double xhas, double yhas, double beta,bool relativitisk_newton){
  bodies.push_back(CelestialBody(masse, xpos, ypos,xhas,yhas,beta,relativitisk_newton));
  return bodies.back(); //returnerer siste element i vektoren som er nylig addet
}

int solsystem::numberOfBodies(){
	return bodies.size();
}
void solsystem::printBodies() {
	for (int i = 0; i < bodies.size(); i++) {
		cout << "Body " << i << ": " << bodies[i].masse << endl;
	}
}
void solsystem::kalkulering_akselerasjon(vector<CelestialBody>& bodies){
  double kin,pot,total,ang_m,t = 0.0;
  double FinalTime,Numberofhs; FinalTime = 200; Numberofhs = 10000;
  double h = FinalTime/Numberofhs;
  ofstream fil;
  fil.open("Energi.txt", ios::app);
  for (int i = 0;i <numberOfBodies(); i++){
    CelestialBody body1 = bodies[i];
    for(int j = 0; j <numberOfBodies();j++ ){
      if (j != i){
        CelestialBody body2 = bodies[j];
        double r12 = sqrt(pow((body1.xpos - body2.xpos), 2) + pow((body1.ypos - body2.ypos), 2));
        body1.xaks = body1.xaks - FourPi2*body2.masse*(body1.xpos - body2.xpos) / pow(r12, body1.beta+1);
        body2.xaks = body2.xaks + FourPi2*body2.masse*(body1.xpos - body2.xpos) / pow(r12, body1.beta+1);
        body1.yaks = body1.yaks - FourPi2*body2.masse*(body1.ypos - body2.ypos) / pow(r12, body1.beta+1);
        body2.yaks = body2.yaks + FourPi2*body2.masse*(body1.ypos - body2.ypos) / pow(r12, body1.beta+1);
        bodies[j] = body2;
      }
    }
    bodies[i] = body1;
    pot = -FourPi2 * bodies[i].masse / bodies[i].r;
    kin = 0.5*bodies[i].masse * (bodies[i].xhas*bodies[i].xhas + bodies[i].yhas*bodies[i].yhas);
    total = pot + kin;
    ang_m = sqrt(pow(bodies[i].xpos*bodies[i].yhas, 2) + pow(bodies[i].ypos*bodies[i].xhas, 2));
    t ++;
    fil << setiosflags(ios::showpoint | ios::uppercase);
    fil << setw(20) << setprecision(8) << t ;
    fil << setw(20) << setprecision(8) << kin;
    fil << setw(20) << setprecision(8) << pot;
    fil << setw(20) << setprecision(8) << total;
    fil << setw(20) << setprecision(8) << ang_m
    << endl;
  }
}
void solsystem::kjoring_algoritme(CelestialBody& body, string outfilename, double FinalTime, int Numberofhs,bool valg_av_algortime,bool relativitisk_newton){
  // ofstream toFile;
  // toFile.open(outfilename, ios::app);
  kalkulering_akselerasjon(bodies);
	if (valg_av_algortime){
    euler(body, FinalTime, Numberofhs,outfilename);
	}
	else{
    //velg 0 for generel og 1 for relativitisk_newton
		velocityVerlet(body, FinalTime, Numberofhs,outfilename,relativitisk_newton);
	}
	// toFile << setprecision(8) << fixed << body.xpos << setw(20) << body.ypos << endl;
}
void solsystem::merkur_presesjon (CelestialBody& body, double FinalTime,int Numberofhs,string outfilename){
  ofstream merkur_fil;
  merkur_fil.open("Merkur.txt");
  ofstream presesjon_fil;
  presesjon_fil.open("presesjon.txt");
  int counter = 0;
  for(int i = 1; i < Numberofhs;i++ ){
    kalkulering_akselerasjon(bodies);
    velocityVerlet(body, FinalTime, Numberofhs,outfilename,1);
    merkur_fil << setprecision(8) << fixed << body.xpos << setw(20) << body.ypos << endl;
    //(body.r*body.r* - 0.3075*0.3075)  < 1E-12
    if(body.r < 0.3075){
      presesjon_fil << setprecision(8) << fixed <<counter <<  setw(20) << atan(body.ypos / body.xpos) <<endl;
      counter ++;
    }
  }
}


void solsystem::beregning_tid(solsystem& system) {

	int FinalTime = 20;
	ofstream myFile;
	myFile.open("Tid.txt");
	myFile << setprecision(8) << fixed << "N " << setw(20) << "Euler(s) " << setw(20) << "Verlet(s)" << endl;

	for (int i = 10e5; i > 0 ; i=i/10) {

		high_resolution_clock::time_point te1 = high_resolution_clock::now();
		system.kjoring_algoritme(system.bodies[0], "Earth.txt", FinalTime, i, 1,0);
		high_resolution_clock::time_point te2 = high_resolution_clock::now();
		duration<double> eulerTime = duration_cast<duration<double>>(te2 - te1);

		high_resolution_clock::time_point tv1 = high_resolution_clock::now();
		system.kjoring_algoritme(system.bodies[0], "Earth.txt", FinalTime, i, 0,0);
		high_resolution_clock::time_point tv2 = high_resolution_clock::now();
		duration<double> verletTime = duration_cast<duration<double>>(tv2 - tv1);

		myFile << setprecision(8) << fixed << i << setw(20) << eulerTime.count() << setw(20) << verletTime.count() << endl;
	}
}
