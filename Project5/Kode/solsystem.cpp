#include"solsystem.h"
#include"celestialBody.h"
#include "algoritmer.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <math.h>
#include<vector>
#include<string>

using namespace std;

solsystem::solsystem(){

}

CelestialBody& solsystem::lag_body(double masse, double xpos, double ypos, double xhas, double yhas, double beta){
  bodies.push_back(CelestialBody(masse, xpos, ypos,xhas,yhas,beta));
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
  for (int i = 0;i <numberOfBodies(); i++){
    CelestialBody body1 = bodies[i];
    for(int j = 0; j <numberOfBodies();j++ ){
      if (j != i){
        CelestialBody body2 = bodies[j];
        double r12 = sqrt(pow((body1.xpos - body2.xpos), 2) + pow((body1.ypos - body2.ypos), 2));
        body1.xaks = body1.xaks - FourPi2*body2.masse*(body1.xpos - body2.xpos) / pow(r12, body1.beta+1);
        body2.xaks = body2.xaks + FourPi2*body2.masse*(body1.xpos - body2.xpos) / pow(r12, body1.beta+1);

        body1.yaks = body1.yaks - FourPi2*body2.masse*(body1.xpos - body2.xpos) / pow(r12, body1.beta+1);
        body2.yaks = body2.yaks + FourPi2*body2.masse*(body1.xpos - body2.xpos) / pow(r12, body1.beta+1);
        bodies[j] = body2;
      }
    }
    bodies[i] = body1;
  }
}
void solsystem::kjoring_algoritme(CelestialBody& body, string outfilename, double FinalTime, int Numberofhs,bool valg_av_algortime){
  kalkulering_akselerasjon(bodies);
	if (valg_av_algortime){
    euler(body, FinalTime, Numberofhs,outfilename );
	}
	else{
		velocityVerlet(body, FinalTime, Numberofhs,outfilename);
	}

}
