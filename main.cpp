#include <iostream>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <string>
#include "time.h"
#include <fstream>



using namespace std;
ofstream ofile;
inline double f(double x) { return 100.0 * exp(-10.0 * x); };
void gauss_generel(int n, double* a, double* b, double* c, double* u, double* s);
void gauss_special(int n, double* b, double* u, double* s);

int main(int argc, char* argv[]) {
	int ex; int test;
	string filename;

	if (argc <= 1) {
		cout << "Bad Usage: you forgot to enter values."  << endl;
		exit(1);
	}
	else {
		ex = atoi (argv[1]);
		test = atoi(argv[2]);
		filename = argv[3];
	}


	for (int i = 1; i <= ex; i++) {
		int  n = (int)pow(10.0, i);
		string fileout = filename;
		string argument = to_string(i);
		fileout.append(argument);


		double h = 1.0 / n; //stepzise
		int n1 = n + 1; int n2 = n + 2;
		double* a = new double[n1];
		double* b = new double[n2];
		double* c = new double[n1];
		double* u = new double[n2];
		double* s = new double[n2];

		//defining the values of elemens
		for (int i = 0; i < n + 2; i++) {
			a[i] = -1.0;
			b[i] = 2.0;
			c[i] = -1.0;
			s[i] = h*h*f(h * i);
		}
		u[0] = u[n2] = 0.0; //boundary condition
		b[n2] = 2.0;


		//couting the time
		clock_t start, finish;
		if (test == 0) {
			start = clock();
			gauss_generel(n, a, b, c, u, s);
			finish = clock();
			double timeSpent = (double(finish - start) / CLOCKS_PER_SEC);
			cout << "Time for generel: " << scientific << timeSpent << setprecision(10) << " s med n = " << n<< endl;

		}
		else if (test == 1) {
			start = clock();
			gauss_special(n, b, u, s);
			finish = clock();
			double timeSpent = (double(finish - start) / CLOCKS_PER_SEC);
			cout << "Time for special: " << scientific << timeSpent << setprecision(10) << " s med n = " << n<< endl;
		}

		/*//writing to files
		ofile.open(fileout);
		ofile << setiosflags(ios::showpoint | ios::uppercase);
		for (int i = 1; i < n; i++) {
			ofile << setw(15) << setprecision(8) << u[i] << endl;
		}
		ofile.close();
		delete[] u;*/
	}


	return 0;
}

void gauss_generel(int n, double* a, double* b, double* c, double* u, double* s)
{
	for (int i = 2; i < n + 1; i++) {
		b[i] = b[i] - ((a[i - 1] * c[i - 1]) / b[i - 1]);
		s[i] = s[i] - ((s[i - 1] * a[i - 1]) / b[i - 1]);
	}
	u[n + 1] = s[n + 1] / b[n + 1];
	for (int i = n + 1; i > 1; i--)
	{
		//u[i] = (s[i] - c[i] * u[i + 1]) / b[i];
		u[i - 1] = (s[i - 1] - c[i - 1] * u[i]) / b[i - 1];
	}
}
void gauss_special(int n, double* b, double* u, double* s)
{
	for (int i = 2; i < n + 1; i++) {
		b[i] = (i + 1.0) / (i);
		s[i] = s[i] + (s[i - 1] / b[i - 1]);
	}
	for (int i = n + 1; i > 1; i--)
	{
		//u[i] = (s[i] + u[i + 1])/ b[i];
		u[i - 1] = (s[i - 1] + u[i]) / b[i - 1];
	}

}
