#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>


using namespace std;
inline double f(double x) { return 100.0 * exp(-10.0 * x); }
void gauss_generel(int n, double* a, double* b, double* c, double* u, double* s);
void gauss_special(int n, double* b, double* u, double* s);

int main() {
    int n;
    cin >> n;
    int n1 = n + 1; int n2 = n + 2;
    double* a = new double[n1];
    double* b = new double[n2];
    double* c = new double[n1];
    double* u = new double[n2];
    double* s = new double[n2];
    double h = 1 / n; // stepzise h = (xn - x0)/n

    //defining the values of elemens
    for (int i = 0; i < n + 1; i++) {
        a[i] = - 1.0;
        b[i] = 2.0;
        c[i] = - 1.0;
        s[i] = h * h * f(h*i);
    }
    u[0] = u[n2] = 0; //boundary condition
    b[n2] = 2.0;


    return 0;
}

void gauss_generel(int n, double* a, double* b, double* c, double* u, double* s)
{
    for (int i = 2; i < n +1; i++) {
        b[i] = b[i] - ((a[i - 1] * c[i - 1]) / b[i - 1]);
        s[i] = s[i] - ((s[i - 1] * a[i - 1]) / b[i - 1]);
    }
    u[n - 1] = s[n - 1] / b[n - 1];
    for (int i = n-2; i < 1; i--)
    {
        u[i] = (s[i] - c[i] * u[i + 1]) / b[i];
    }
}
void gauss_special(int n, double* b, double* u, double* s)
{
    for (int i = 1; i <= n+1; i++){
        b[i] = (i + 1.0)/(i);
        s[i] = s[i] + (s[i-1]/b[i-1]);
    }
    for (int i = n-2; i < 1; i--)
    {
        u[i] = (s[i] + u[i + 1])/ b[i];
    }

}
/*
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;


inline double f(double x){return 100.0*exp(-10.0*x);} //function we want to solve  --RHS

void general_solver(int, double*, double*, double*, double*, double*);

int main(int argc, char* argv[]){
    int N;
    cin >> N;

    double dt = 1.0/(N+1); //time step

    //allocating memory dynamicly for our vectors which will be used to calculate solution
    double* a = new double[N+1]; //sub main diagonal
    double* b = new double[N+2]; //main diagonal
    double* c = new double[N+1]; //above main diagonal
    double* s = new double[N+2]; //our RHS function --f(x)*h^2
    double* v = new double[N+2]; //simulated solution


    //initializing the arrays/vectors
    for(int i = 0; i < N+2; i++){
        s[i] = dt*dt*f(dt*i);
        a[i] = -1.0;
        b[i] = 2.0;
        c[i] = -1.0;
    }
    //Dirichlet boundary conditions
    s[0]= s[N+2] = 0;

    //the last main diagonal element
    b[N+2] = 2.0;

    for(int i =1; i < N+1; i++){
        cout << v[i] << endl;
    }

}

//general gaussian elimination - all the elements can have different different values
void general_solver(int N, double* a, double* b, double* c, double* s, double* v){
    //forward
    for(int i = 2; i <= N+1; i++){
        b[i] -= (a[i-1]*c[i-1])/b[i-1];
        s[i] -= s[i-1]*(a[i-1]/b[i-1]);

    }
    //backward
    v[N+1] = s[N+1]/b[N+1]; //initial value (last element)

    for(int i = N+1; i > 1; i--){
        v[i-1]= (s[i-1] - c[i-1]*v[i])/b[i-1];

    }
}

*/
