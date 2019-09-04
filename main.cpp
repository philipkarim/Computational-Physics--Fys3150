#include <iostream>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <string>
#include "time.h"

using namespace std;
inline double f(double x) { return 100.0 * exp(-10.0 * x); };
void gauss_generel(int n, double* a, double* b, double* c, double* u, double* s);
void gauss_special(int n, double* b, double* u, double* s);

int main() {
    int n; int test;
    cin >> n;
    cin >> test;
    int n1 = n + 1; int n2 = n + 2;
    double* a = new double[n1];
    double* b = new double[n2];
    double* c = new double[n1];
    double* u = new double[n2];
    double* s = new double[n2];
    double h = 1.0 / n; // stepzise h = (xn - x0)/n

    //defining the values of elemens
    for (int i = 0; i < n + 2; i++) {
        a[i] = - 1.0;
        b[i] = 2.0;
        c[i] = - 1.0;
        s[i] = h * h * f(h*i);
    }
    u[0] = u[n2] = 0.0; //boundary condition
    b[n2] = 2.0;

    //couting the time
    clock_t start, finish;
    if(test == 0){
        start = clock();
        gauss_generel(n,a,b,c,u,s);
        finish = clock();
        double timeSpent = (double(finish - start)/CLOCKS_PER_SEC );
        cout << "Time spent for generel: " << fixed << timeSpent <<setprecision(5)<< " s."<< endl;
    } else if (test == 1) {
        start = clock();
        gauss_special(n,b,u,s);
        finish = clock();
        double timeSpent = (double(finish - start)/CLOCKS_PER_SEC );
        cout << "Time spent for special: " << fixed << timeSpent <<setprecision(5)<< " s."<< endl;
    }else {
        cout << "Velg 0 for den generelle Gauss eller 1 for den spesielle Gauss." << endl;
        exit(1)
    }

    return 0;
}

void gauss_generel(int n, double* a, double* b, double* c, double* u, double* s)
{
    for (int i = 2; i < n +1; i++) {
        b[i] = b[i] - ((a[i - 1] * c[i - 1]) / b[i - 1]);
        s[i] = s[i] - ((s[i - 1] * a[i - 1]) / b[i - 1]);
    }
    u[n + 1] = s[n + 1] / b[n + 1];
    for (int i = n+1; i > 1; i--)
    {
        //u[i] = (s[i] - c[i] * u[i + 1]) / b[i];
        u[i-1]= (s[i-1] - c[i-1]*u[i])/b[i-1];
    }
}
void gauss_special(int n, double* b, double* u, double* s)
{
    for (int i = 2; i < n+1; i++){
        b[i] = (i + 1.0)/(i);
        s[i] = s[i] + (s[i-1]/b[i-1]);
    }
    for (int i = n+1; i > 1; i--)
    {
        //u[i] = (s[i] + u[i + 1])/ b[i];
        u[i-1]= (s[i-1] + u[i])/b[i-1];
    }

}
