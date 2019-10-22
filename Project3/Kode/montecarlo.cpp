#include <iostream>
#include <cmath>
#include <omp.h>
#include <armadillo>


//importerer disse for å beregne tiden:
#include <ctime>
#include <ratio>
#include <chrono>

using namespace std::chrono; //beregner tiden med chrono
using namespace std;
using namespace arma;

double integrate(double x1, double y1, double z1, double x2, double y2, double z2);
double ingerate_impor(double r1, double r2, double costheta1, double costheta2, double phi1, double phi2);
double relativError(double nummerical, double exact);

int main() {
	long int n; double a, b;
	n = 10e5; a = -3.0; b = 3.0; //bestemmer N og grensene her

	double fx, sum_sigma, crude_mc, var_mc, int_mc,int_exact,rela_feil_mc, rela_feil_mc_imSa;
	double r1,r2,costheta1,costheta2,phi1,phi2,l, var_mc_imSa, int_mc_imSa;
	const double PI = 3.141592653589793238463;

	vec y = zeros<vec>(6);
	int_exact = (5 * PI * PI) / (16 * 16);

	//BRUTE FORCE
	arma_rng::set_seed_random(); // forandrer RNG seed
	crude_mc = 0.0; fx = 0.0; sum_sigma = 0.0;
	high_resolution_clock::time_point t1 = high_resolution_clock::now(); //Start tid
	for (int i = 1; i <= n; i++) {
		//fyller x med 6 random verdier fra 0 til 1.
		vec x = randu<vec>(6);

		//beregner y verdiene
		y = a + x * (b - a);

		fx = integrate(y(0), y(1), y(2), y(3), y(4), y(5));
		crude_mc += fx;
		sum_sigma += fx * fx;
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();//slutt tid

	crude_mc = crude_mc/((double) n);
	sum_sigma = sum_sigma/((double)n);
	var_mc = sum_sigma - crude_mc* crude_mc;
	int_mc = crude_mc * pow((b - a), 6);

	//beregner tidsforskjellen
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

	//beregner relativ feil
	rela_feil_mc = relativError(int_mc, int_exact);
	cout << "N  " << n<< endl;

	cout << "Brute force  " << endl;
	cout << "I: " << int_mc <<endl;
	cout << "Exact: " << int_exact << endl;
	cout << "Relativ feil: " << rela_feil_mc << endl;
	cout << "Standard deviation: " << pow((b - a), 6) * sqrt(var_mc / ((double)n)) <<endl;
	cout << "Tid: " << time_span.count() << " s" << endl;

	//IMPORTANCE SAMPLING
	arma_rng::set_seed_random(); // forandrer RNG seed
	high_resolution_clock::time_point t3 = high_resolution_clock::now();//start tid
	crude_mc = 0.0; fx = 0.0; sum_sigma = 0.0; l = 1.0 / 4;
	//paralellisering:
	//#pragma omp parallel for private(i) reduction(+:crude_mc) reduction(+:sum_sigma)
	for (int i = 1; i <= n; i++) {
		vec x1 = randu<vec>(6);

		r1 = -l * log(1 - x1(0));
		r2 = -l * log(1 - x1(1));
		costheta1 = 2 * x1(2) - 1;
		costheta2 = 2 * x1(3) - 1;
		phi1 = 2 * PI * x1(4);
		phi2 = 2 * PI * x1(5);

		fx = ingerate_impor(r1, r2, costheta1, costheta2, phi1, phi2);
		crude_mc += fx;
		sum_sigma += fx * fx;
	}
	high_resolution_clock::time_point t4 = high_resolution_clock::now();//slutt tid

	crude_mc = crude_mc / ((double)n);
	sum_sigma = sum_sigma / ((double)n);
	var_mc_imSa = sum_sigma - crude_mc * crude_mc;
	int_mc_imSa = crude_mc*2*2*(2*PI)*(2 * PI)*l*l;

	rela_feil_mc_imSa = relativError(int_mc_imSa, int_exact);
	duration<double> time_span1 = duration_cast<duration<double>>(t4 - t3);

	cout << "Importance sampling " << endl;
	cout << "I: " << int_mc_imSa << endl;
	cout << "Exact: " << int_exact << endl;
	cout << "Relativ feil: " << rela_feil_mc_imSa << endl;
	cout << "Standard deviation: " << 2*2*(2 * PI) * (2 * PI) * l * l* sqrt(var_mc_imSa /n) << endl;
	cout << "Tid: " << time_span1.count() << " s" << endl;

	return 0;
}

double integrate(double x1, double y1, double z1, double x2, double y2, double z2){
	int alpha = 2;
	double r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
	double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);

	double ledd1= exp(-2 * alpha * (r1 + r2));
	double ledd2 = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2));

	if (ledd2 < 1e-10){ //unngå å dele på null
		return 0;
	}
	else {
		return (ledd1 / ledd2);
	}
}
double ingerate_impor(double r1,double r2,double costheta1,double costheta2,double phi1,double phi2){
	double ledd1 = r1 * r1 * r2 * r2;
	double cosB = costheta1 * costheta2 + sqrt(1 - costheta1 * costheta1) * sqrt(1 - costheta2 * costheta2) * cos(phi1 - phi2);
	double ledd2 = sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * cosB);
	if (ledd2 < 1e-10) { //unngå å dele på null
		return 0;
	}
	else {
		return (ledd1 / ledd2);
	}
}
double relativError(double nummerical, double exact) {
	double r;
	double epslon = 1e-10;
	if (exact <= epslon) { //lite test for å ungå å dele på null.
		r = 0;
	}
	else {
		r = abs((nummerical - exact) / (exact)); //relativ feil.
	}
	return r;
}
