/*#include <iostream>
using namespace std;

int main()
{
    double a, b;

    a = 3.6;
    b = 2.5;
    cout << "Hello World" << endl;
    cout << a + b;
    return 0;
}
*/
/* Test for å fikse armadillo
#include <iostream>
#include <armadillo>

using namespace std;

int main()
{
    cout << "Hello World!" << endl;

    unsigned int n = 10;
    arma::Mat<double> A(n, n, arma::fill::zeros);

    for(unsigned int i=0; i<n; i++){
        A(i,i) = 10;
    }

    for(unsigned int i=0; i<n-3; i++){
        A(i+3,i) = -5;
        A(i,i+3) = -5;
    }

    arma::Mat<double> P, L, U;

    arma::lu(P, L, U, A);

    cout << A << endl << P << endl << L << endl << U << endl;

    return 0;
}*/

//Importerer nødvendige pakker
#include <iostream>
#include <cmath>

using namespace std;
//Denne kommer vi til å måtte bruke i d(?)
//using namespace arma;

//Funksjonen og den dobbeltderiverte(Til bruk i error kalkulasjoner)
inline double funksjon(double x){
    return 100.0*exp(-10.0*x);
}
inline double derivert2(double x) {
    return 1.0-(1-exp(-10))*x-exp(-10*x);
}

//Main funksjonen
int main() {
    int n = 5;

//Definerer vektorer. Usikker på om dette er riktig definert siden pcen
//min tydeligvis hater alt som har med armadillo å gjøre
    double *d = new double[n];
    double *d_tilde = new double[n];
    double *g = new double[n];
    double *g_tilde = new double[n];
    double *a = new double[n];
    double *b = new double[n];
    double *u = new double[n];
    double *q = new double[n];

    return 0;
}


//Her vettu, her er forward og bacward substitution
//Litt usikker på om substitutionene dine var riktige så tok den fra kompendiet
void gauss(int n,double* d,double* d_tilde,double* g,double* g_tilde,double* a,double*b,double* u,double*q)
{
    //Den første forward i koden
    d[1] = d_tilde[1];
    for (int i = 2; i < n - 1; i++) {

        d_tilde[i] = d[i] - (a[i - 1] * b[i - 1]) / d_tilde[i - 1];
        g_tilde[i] = g[i] - (g_tilde[i - 1] * a[i - 1]) / d_tilde[i - 1];
    }
    //Forward fra kompendiet(s186) med våre verdier (måtte legge til en array t
    //d er diagonal, a er diagonal over d, b er diagonal under d, u er v vektroen, g er svaret
    double qt = d[1];
    u[1] = g[1]/qt;
    for(int i=2 ; i <= n ; i++) {
        q[i] = a[i-1]/qt;
        qt = d[i]-b[i]*q[i];
        u[i] = (g[i] - b[i]*u[i-1])/qt;
    }

    //Den første backward i koden
    u[n - 1] = g_tilde[n - 1] / d_tilde[n - 1];
    for (int i = 0; i < n - 1; i++)
    {
        u[i] = (g_tilde[i] - b[i] * u[i + 1]) / d_tilde[i];
    }

    //Backward fra kompendiet med våre verdier og vektorer
    for(int i=n-1 ; i >= 1 ; i--) {
    u[i] -= q[i+1]*u[i+1];
    }



}
