#include "algoritmer.h"
#include "celestialBody.h"
#include "solsystem.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include<string>

using namespace std;
double FourPi2 = 4*pi*pi;


void euler(CelestialBody& body, double FinalTime, int Numberofhs,string outfilename){
  double h = FinalTime/Numberofhs;
  double kinetisk_energi,potensiell_energi, total_energi,angular_moment;
  double time =0;
  //loser diffligningen:
  body.r = sqrt(body.xpos*body.xpos + body.ypos*body.ypos);
  body.xpos += h*body.xhas;
  body.ypos += h*body.yhas;
  body.xhas -= h*FourPi2*body.xpos/(body.r*body.r*body.r);
  body.yhas -= h*FourPi2*body.ypos/(body.r*body.r*body.r);
  body.r = sqrt(body.xpos*body.xpos+body.ypos*body.ypos);
  //energi:
  kinetisk_energi = FourPi2 *0.5*body.masse*(body.xhas +body.yhas)*(body.xhas+body.yhas);
  potensiell_energi = -FourPi2*FourPi2*body.masse/(body.r);
  total_energi = kinetisk_energi+potensiell_energi;
  angular_moment += sqrt(pow(body.xpos*body.yhas, 2) + pow(body.ypos*body.xhas, 2));
  time += h;
  //skriver til fil
  output(time, body.xpos, body.ypos, body.xhas, body.yhas,outfilename,kinetisk_energi,potensiell_energi,total_energi,angular_moment);   // write to file
}

void velocityVerlet(CelestialBody& body, double FinalTime, int Numberofhs,string outfilename,bool relativitisk_newton) {
    double kinetisk_energi,potensiell_energi, total_energi,angular_moment;
    double nest_xaks,nest_yaks;
    double h = FinalTime/Numberofhs;
    double time = 0;
    //loser diffligningen:
    body.xpos = body.xpos+h*body.xhas +(h*h/2)*body.xaks;
    body.ypos =  body.ypos+h*body.yhas+(h*h/2)*body.yaks;
    body.r = sqrt(body.xpos*body.xpos+body.ypos*body.ypos);
    if(!relativitisk_newton){
      nest_xaks = -(FourPi2*body.xpos) / (pow(body.r, body.beta + 1));
      nest_yaks = -(FourPi2*body.ypos) / (pow(body.r, body.beta + 1));
    }else{
        body.xaks = -((FourPi2*body.xpos) / (pow(body.r, body.beta + 1))) *(1+3*(pow(body.xpos*body.yhas,2) + pow(body.ypos*body.xhas,2)))/(pow(body.r,2)*pow(3e8*60*60*24*365/ 1.496e8, 2));
        body.yaks = -((FourPi2*body.ypos) / (pow(body.r, body.beta + 1))) *(1+3*(pow(body.xpos*body.yhas,2) + pow(body.ypos*body.xhas,2)))/(pow(body.r,2)*pow(3e8*60*60*24*365/ 1.496e8, 2));
    }

    body.xhas= body.xhas+(h/2)*body.xaks+(h/2)*nest_xaks;
    body.yhas= body.yhas+(h/2)*body.yaks+(h/2)*nest_yaks;
    body.xaks= nest_xaks;
    body.yaks=nest_yaks;
    time +=h;
    //energi:
    kinetisk_energi += FourPi2 *0.5*body.masse*(body.xhas +body.yhas)*(body.xhas+body.yhas);
    potensiell_energi = - FourPi2*FourPi2*body.masse/(body.r);
    total_energi = kinetisk_energi+potensiell_energi;
    angular_moment += sqrt(pow(body.xpos*body.yhas, 2) + pow(body.ypos*body.xhas, 2));
    //skriver til fil
    output(time, body.xpos, body.ypos, body.xhas, body.yhas,outfilename,kinetisk_energi,potensiell_energi,total_energi,angular_moment);

}
void output(double time, double x, double y, double xhas, double yhas,string outfilename,double kinetisk_energi,double potensiell_energi,
  double total_energi,double angular_moment)
{
  ofstream ofile;
  ofile.open(outfilename, ios::app);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(20) << setprecision(8) << time ;
  ofile << setw(20) << setprecision(8) << x;
  ofile << setw(20) << setprecision(8) << y;
  ofile << setw(20) << setprecision(8) << xhas;
  ofile << setw(20) << setprecision(8) << yhas;
  ofile << setw(20) << setprecision(8) << kinetisk_energi;
  ofile << setw(20) << setprecision(8) << potensiell_energi;
  ofile << setw(20) << setprecision(8) << total_energi;
  ofile << setw(20) << setprecision(8) << angular_moment
  << endl;
}
