#include <iostream>
#include <armadillo>
#include <iomanip>
#include <cmath>
#include <cstdlib>
//#include "jacobimetode.h"
#include "prosjekt2funksjoner.h"
using namespace std;
using namespace arma;

int n=3;
double rhomin=0.0;
double rhomax=5.0;
double h = (rhomax-rhomin)/(n);

double d = 2.0/(h*h);
double a = -1.0/(h*h);

vec eigval_arma;
vec eigval_jacobi;

int z;

//Hovedmatrisen vår A, test matrisen vår B (random inputs), Identitetsmatrisen vår C som skal fylles med egenverdier

    mat A =zeros<mat>(n,n);//Jacobi
    mat E =zeros<mat>(n,n);//Schrodinger
    mat B = mat(n,n,fill::randu);
    mat D=randu<mat>(n,n);

    mat C= eye<mat>(n,n);

    //vec e er første kolonnen, kan skrive toeplitz(e,e2) for kolonne og rad
    //kan antagelig bruke vektorene fra prosjekt 1 til å lage matrisen vår
    vec e = randu<vec>(n);
    mat F = toeplitz(e);


//laget en funksjon som fyller matrisen vår
//denne er kopiert
void makematrix1(int n, double d,double a){
    for(int i = 0; i < n; i++){
        A[i,i] = d;
        if(i != n-1 ){
            A[i+1,i] = a;
            A[i,i+1] = a;
            }}}

//må lage en funksjon her for å fylle matrise når det er potensial tilgjengelig
//void makematrix2




int main()
{
   cin>>z>>endl;

   if (z==1){
        //Armadillo
        mat C = diagmat(A);     //Diagonaliserer
        //cout<<B<<endl;
        //cout<<D<<endl;
        //Finner egenvektorer og egenverdier
        eig_sym(eigval_arma, C);

        cout<<eigval_arma<<endl;

    if(z==2){
        jacobi_method(A, C, n);
    }

    }

    return 0;
}
