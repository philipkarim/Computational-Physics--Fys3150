#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
//importerer disse for å beregne tiden:
#include <ctime>
#include <ratio>
#include <chrono>

using namespace std;
ofstream ofile;
//definerer funksjonene:
inline double f(double x) { return 100.0 * exp(-10.0 * x); };
inline double exact(double x) { return (1.0 - (1 - exp(-10)) * x - exp(-10 * x)); }

void gauss_generel(int n, double* a, double* b, double* c, double* u, double* s);
void gauss_special(int n, double* b, double* u, double* s);
double relativError(double nummerical, double exact);

int main(int argc, char* argv[]) {
	/*brukeren må oppgi tre argumenter. forste er eksponenent, andre: 0 (generell) eller 1 (spesial), tredje: navn på filen.
	*/
	int ex; int test;
	string filename;
	
	if (argc <= 1) {
		cout << "Bad Usage: you forgot to enter values." << endl;
		exit(1);
	}
	else {
		ex = atoi(argv[1]); //eksponent
		test = atoi(argv[2]); //0 eller 1
		filename = argv[3];//filnavn
	}

	for (int i = 1; i <= ex; i++) {
		int  n = (int)pow(10.0, i); //regner ut n, avhengig av hva brukeren har skrevet.
		//henter navnet til filen, og lager n filer.
		string fileout = filename;
		string argument = to_string(i);
		fileout.append(argument);
		
		double h = 1.0 / (n+1); //stepzise
		int n1 = n + 1; int n2 = n + 2;
		//definerer arrays
		double* a = new double[n1];
		double* b = new double[n2];
		double* c = new double[n1];
		double* u = new double[n2];
		double* s = new double[n2];
		double* error = new double[n2];


		//setter inn verdiene til arrays
		for (int i = 0; i < n + 2; i++) {
			a[i] = -1.0;
			b[i] = 2.0;
			c[i] = -1.0;
			s[i] = h * h * f(h * i);
		}
		u[0] = u[n1] = 0.0; //grensebetingelsen.
		b[n2] = 2.0; //siste plassen i b 

		//beregner tiden med chrono
		using namespace std::chrono;
		//bestemmer om å regne generelle eller spesielle
		if (test == 0) {
			high_resolution_clock::time_point t1 = high_resolution_clock::now(); //Start tid
			gauss_generel(n, a, b, c, u, s);
			high_resolution_clock::time_point t2 = high_resolution_clock::now();//slutt tid
			//beregner tidsforskjellen
			duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

			cout << "Time for generel: " << time_span.count() << " s med n = " << n << endl;

		}
		else if (test == 1) {
			high_resolution_clock::time_point t1 = high_resolution_clock::now(); //Start tid
			gauss_special(n, b, u, s);
			high_resolution_clock::time_point t2 = high_resolution_clock::now();//slutt tid
			//beregner tidsforskjellen
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
			//fyller error med ulike verdier for relative feil
			error[i] = relativError(u[i], exact(h*i)); 
			//beregner storste feilen:
			if (error[i] >= abs(max_error)) {
				max_error = error[i];
			}
		}
		
		//skriver til andre filer
		ofile.open(fileout);
		ofile << setiosflags(ios::showpoint | ios::uppercase);
		//sender error: log10(h) --- log10(error)
		/*ofile << setprecision(12) << left << setw(25) << log10(h)
			<< left << setw(25) <<  log10(max_error)
			<< endl;*/
		for (int z = 1; z < n; z++) {
			//sender losningen til fil
			ofile << setw(15) << setprecision(8) << u[z] <<endl;
		}
		ofile.close();
		delete[] u; 
	}
	return 0;
}

void gauss_generel(int n, double* a, double* b, double* c, double* u, double* s)
{	//forward subsitution
	for (int i = 2; i < n + 1; i++) {
		b[i] = b[i] - ((a[i - 1] * c[i - 1]) / b[i - 1]);
		s[i] = s[i] - ((s[i - 1] * a[i - 1]) / b[i - 1]);
	}
	//backward substitution
	for (int i = n+1; i > 1; i--)
	{
		//u[i] = (s[i] - c[i] * u[i + 1]) / b[i];
		u[i - 1] = (s[i - 1] - c[i - 1] * u[i]) / b[i - 1];
	}
}
void gauss_special(int n, double* b, double* u, double* s)
{	//forward subsitution
	for (int i = 2; i < n + 1; i++) {
		b[i] = (((double)i) + 1.0) / ((double)i);
		s[i] = s[i] + (s[i - 1] / b[i - 1]);
	}
	//backward substitution
	u[n + 1] = s[n + 1] / b[n + 1];
	for (int i = n + 1; i > 1; i--)
	{
		u[i - 1] = (s[i - 1] + u[i]) / b[i - 1];
	}

}
double relativError(double nummerical, double exact){
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
