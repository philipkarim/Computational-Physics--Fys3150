#include <iostream>
#include <armadillo>
#include <vector>
#include <array>
#include "prosjekt2funksjoner.h"
#include "unixtest.h"

using namespace std;
using namespace arma;

int main(){

//using namespace std::vector;

//skriver 2 tester, teste at vi får maxdiagonal og teste at vi får riktige egenverdier

//maxdiagonal test

    double epsilon=10e-8;
    int n=5;

    int k,l;
    mat Atest = zeros<mat>(n, n);
    mat Ctest = eye<mat>(n, n);


    int * ar = new int[n];
    ar[0] = 1;
    ar[1] = 2;
    ar[2] = 3;
    ar[3] = 4;
    ar[4] = 5;

    for(int i = 0; i < n; i++ ){
        for (int j = 0; j < n; j++){
            Atest(i,j) = ar[i];
        }
    }
    //setter max verdien vår til 9
    Atest(3,2)=9;
    cout <<Atest;


    maxoffdiag(Atest,k,l,n);
    int max=9;
    if (abs(max-9.0<epsilon)){
        cout<< "Test bestatt";}
    else {
        cout<<"Test ikke bestått";
    }

//Riktige egenverdier test

    double * ev_jacobi_test = new double[n];
    double * ev_jacobi_test2 = new double[n];

    double * ev_fasit_test=new double[n];

    jacobi_method (Atest,Ctest,n);

    for(int i = 0; i < n; i++ ){
            Ctest(i,i) = ev_jacobi_test[i];
        }


    //sorterer elementene i stigende rekkefølge
    vec ev_jacobi_test2=sort(ev_jacobi_test);

    //for (int i = 0; i < n; ++i){
      //      cout << ar[i] << ";}

    //selectionsort(ev_jacobi_test,n);

    ev_fasit_test[0]=-0.94097;
    ev_fasit_test[1]=0;
    ev_fasit_test[2]=0;
    ev_fasit_test[3]=0;
    ev_fasit_test[4]=15.940971508067074;

    if (abs(ev_fasit_test[0]-ev_jacobi_test[0])<epsilon &&
        abs(ev_fasit_test[1]-ev_jacobi_test[1])<epsilon &&
        abs(ev_fasit_test[2]-ev_jacobi_test[2])<epsilon &&
        abs(ev_fasit_test[3]-ev_jacobi_test[3])<epsilon &&
        abs(ev_fasit_test[4]-ev_jacobi_test[4])<epsilon){
        cout<<"Test bestått";}
    else {
        cout<<"Test ikke bestått";
    }

}
