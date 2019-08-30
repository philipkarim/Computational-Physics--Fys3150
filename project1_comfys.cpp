#include <iostream>
using namespace std;

int main() {
	int n = 5;

	double *d = new double[n];
	double *d_tilde = new double[n];
	double *g = new double[n];
	double *g_tilde = new double[n];
	double *a = new double[n];
	double *b = new double[n];
	double *u = new double[n];


	return 0;
}

void gauss(int n,double* d,double* d_tilde,double* g,double* g_tilde,double* a,double*b,double* u)
{
	d[1] = d_tilde[1];
	for (int i = 2; i < n - 1; i++) {
		d_tilde[i] = d[i] - (a[i - 1] * b[i - 1]) / d_tilde[i - 1];
		g_tilde[i] = g[i] - (g_tilde[i - 1] * a[i - 1]) / d_tilde[i - 1];
	}

	u[n - 1] = g_tilde[n - 1] / d_tilde[n - 1];
	for (int i = 0; i < n - 1; i++)
	{
		u[i] = (g_tilde[i] - b[i] * u[i + 1]) / d_tilde[i];
	}
}
