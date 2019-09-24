#include "prosjekt2funksjoner.h"

#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

//skriver 2 tester, teste at vi får maxdiagonal og teste at vi får riktige egenverdier

//Definerer vektorer
static vec testeigv_arma;
static vec testeigv_jacobi;

static int n=4;


mat A =zeros<mat>(n,n);
mat B = mat(n,n,fill::randu);
mat D=randu<mat>(n,n);

mat C= eye<mat>(n,n);

    //vec e er første kolonnen, kan skrive toeplitz(e,e2) for kolonne og rad
    //kan antagelig bruke vektorene fra prosjekt 1 til å lage matrisen vår
vec e = randu<vec>(n);
mat F = toeplitz(e);


