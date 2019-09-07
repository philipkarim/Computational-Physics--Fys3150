//LU dekomposisjon og losning av likningen.

#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <ratio>
#include <chrono>
#include <armadillo>
#include <fstream>


using namespace std;
using namespace arma;
using namespace std::chrono;

ofstream ofile;

inline double f(double x) { return 100.0 * exp(-10.0 * x); };
inline double exact(double x) { return (1.0 - (1 - exp(-10)) * x - exp(-10 * x)); }
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
	}

	int n = (int)pow(10.0, ex);
	double h = 1.0 / (n + 1);
	vec y = zeros<vec>(n);
	vec x = zeros<vec>(n);
	vec s = zeros<vec>(n);
	mat A = zeros<mat>(n, n);
	mat L, U;
	for (int i = 0; i < n; i++) {
		A(i, i) = 2;
		if (i - 1 >= 0) {
			A(i, i - 1) = -1;
		}
		if (i + 1 <= n - 1) {
			A(i, i + 1) = -1;
		}
		s(i) = h * h * f(h * i);
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	lu(L, U, A); //dekomperer matrise A til L og U.
	y = solve(L, s);
	x = solve(U, y);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

	cout << "Time for LU: " << time_span.count() << " s med n = " << n << endl;

	vec solution = zeros<vec>(n + 2);
	for (int i = 0; i < n; i++) {
		solution(i + 1) = x(i);
	}
	solution.save("LU_plot", raw_ascii);

	double* error = new double[n + 2];
	error[0] = 0.0;
	error[n+1] = 0.0;
	double max_error = error[0];
	for (int i = 1; i < n+1; i++) {
		error[i] = relativError(solution[i], exact(h * i));
		if (error[i] >= abs(max_error)) {
			max_error = error[i];
		}
	}
	ofstream myfile("LU_error_plot");
	if (myfile.is_open())
	{
		myfile << setprecision(12) << left << setw(25) << log10(h)
			<< left << setw(25) << log10(max_error)
			<< endl;

		myfile.close();
	}

	return 0;
}

double relativError(double nummerical, double exact) {
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
