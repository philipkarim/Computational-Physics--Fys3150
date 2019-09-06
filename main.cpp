#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>

#include <ctime>
#include <ratio>
#include <chrono>

using namespace std;
ofstream ofile;

inline double f(double x) { return 100.0 * exp(-10.0 * x); };
inline double exact(double x) { return (1.0 - (1 - exp(-10)) * x - exp(-10 * x)); }
void gauss_generel(int n, double* a, double* b, double* c, double* u, double* s);
void gauss_special(int n, double* b, double* u, double* s);
double relativError(double nummerical, double exact);

int main(int argc, char* argv[]) {
	int ex; int test;
	string filename;

	if (argc <= 1) {
		cout << "Bad Usage: you forgot to enter values." << endl;
		exit(1);
	}
	else {
		ex = atoi(argv[1]);
		test = atoi(argv[2]);
		filename = argv[3];
	}

	for (int i = 1; i <= ex; i++) {
		int  n = (int)pow(10.0, i);
		string fileout = filename;
		string argument = to_string(i);
		fileout.append(argument);

		double h = 1.0 / (n+1); //stepzise
		int n1 = n + 1; int n2 = n + 2;
		double* a = new double[n1];
		double* b = new double[n2];
		double* c = new double[n1];
		double* u = new double[n2];
		double* s = new double[n2];
		double* error = new double[n2];


		//defining the values of elemens
		for (int i = 0; i < n + 2; i++) {
			a[i] = -1.0;
			b[i] = 2.0;
			c[i] = -1.0;
			s[i] = h * h * f(h * i);
		}
		u[0] = u[n1] = 0.0; //boundary condition
		b[n2] = 2.0;

		//couting the time using chrono

		using namespace std::chrono;
		if (test == 0) {
			high_resolution_clock::time_point t1 = high_resolution_clock::now();
			gauss_generel(n, a, b, c, u, s);
			high_resolution_clock::time_point t2 = high_resolution_clock::now();
			duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

			cout << "Time for generel: " << time_span.count() << " s med n = " << n << endl;

		}
		else if (test == 1) {
			high_resolution_clock::time_point t1 = high_resolution_clock::now();
			gauss_special(n, b, u, s);
			high_resolution_clock::time_point t2 = high_resolution_clock::now();
			duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
			cout << "Time for special: " << time_span.count() << " s med n = " << n << endl;
		}
		else {
			cout << "Error: second argv not valid." << endl;
		}

		//finner den relative feilen:
		error[0] = 0.0;
		error[n1] = 0.0;
		double max_error = error[0];
		for (int i = 1; i < n1; i++){
			error[i] = relativError(u[i], exact(h*i));
			if (error[i] >= abs(max_error)) {
				max_error = error[i];
			}
		}
		/*
		//writing to files
		ofile.open(fileout);
		ofile << setiosflags(ios::showpoint | ios::uppercase);
		//sender error: log10(h) --- log10(error)
		ofile << setprecision(12) << left << setw(25) << log10(h)
			<< left << setw(25) <<  log10(max_error)
			<< endl;
		for (int z = 1; z < n; z++) {
			//ofile << setw(15) << setprecision(8) << u[z] <<endl;
		}
		ofile.close();
		delete[] u; */
	}
	return 0;
}

void gauss_generel(int n, double* a, double* b, double* c, double* u, double* s)
{
	for (int i = 2; i < n + 1; i++) {
		b[i] = b[i] - ((a[i - 1] * c[i - 1]) / b[i - 1]);
		s[i] = s[i] - ((s[i - 1] * a[i - 1]) / b[i - 1]);
	}
	//u[n] = s[n] / b[n];
	for (int i = n+1; i > 1; i--)
	{
		//u[i] = (s[i] - c[i] * u[i + 1]) / b[i];
		u[i - 1] = (s[i - 1] - c[i - 1] * u[i]) / b[i - 1];
	}
}
void gauss_special(int n, double* b, double* u, double* s)
{
	for (int i = 2; i < n + 1; i++) {
		b[i] = (((double)i) + 1.0) / ((double)i);
		s[i] = s[i] + (s[i - 1] / b[i - 1]);
	}
	u[n + 1] = s[n + 1] / b[n + 1];
	for (int i = n + 1; i > 1; i--)
	{
		//u[i] = (s[i] + u[i + 1])/ b[i];
		u[i - 1] = (s[i - 1] + u[i]) / b[i - 1];
	}

}
double relativError(double nummerical, double exact){
	double r;
	double epslon = 1e-10;
	if (exact <= epslon) {
		r = 0;
	}
	else {
		r = abs((nummerical - exact) / (exact));
	}
	return r;
}
