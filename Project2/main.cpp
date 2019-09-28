//JACOBI OPPGAVE 2B and 2D:
#include <iostream>
#include <cmath>
#include <armadillo>
#include <string>
#include <iomanip>
#include <fstream>

using namespace arma;
using namespace std;

ofstream ofile;

void fylleMatriseA(mat& A,int n,double h);
void fylleMatriseA_potensial(mat& A,int n, double h, int rho_0,vec&p);
double maxoffdiag(mat& A, int* k, int* l, int n);
void rotate(mat& A, mat& R, int k, int l, int n);
double* relativError(vec& nummerical, int* exact);

int main(){
	int n = 100; int rho_0 = 0; int rho_N = 10;
	int k, l;
	int z = 1; //z = 0 ingen potensial, z = 1 med potensial 1 e, z = 2 potensial 2 elek
	int w = 0.01;
	double h = (rho_N - rho_0) / ((float)n);
	double epsilon = 1e-10;
	double max_number_iterations = n*n*n;
	int iterations = 0;

	//defining the matrix A and putting in elements:
	mat R = zeros<mat>(n, n);
	mat A = zeros<mat>(n, n);
	vec eig_values = zeros<vec>(n);
	vec eigenvectors = zeros<vec>(n);
	vec arma_eigValue = zeros<vec>(n);
	vec potensial_1e = zeros<vec>(n);
	vec potensial_2e = zeros<vec>(n);
	vec num_eig = zeros<vec>(n);
	int* a_verdi = new int[n];
	a_verdi[0] = 3; a_verdi[1] = 7; a_verdi[2] = 11; a_verdi[3] = 15;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {	//putter inn R elementer
				R(i,j) = 1.0;
			}
			else {
				R(i,j) = 0.0;
			}
		}
		//beregner potensialene for 1 elek og 2 elek:
		potensial_1e(i) = ((rho_0 + (i + 1) * h) * (rho_0 + (i + 1) * h)); //p^2
		potensial_2e(i) = ((rho_0 + (i + 1) * h) * (rho_0 + (i + 1) * h) * w * w) + (1 / (rho_0 + (i + 1) * h)); //w^2*p^2 + 1/p
	}
	//beregner diagonal og ikke-diagonal elementer:
	double d_value = 2.0 / (h * h);
	double ikke_d_value = -1.0 / (h * h);

	//putter inn elementer i A
	if (z == 0){ //ikke potensial
		cout << "Beregner egenverdiene uten potensial:" << endl;
		fylleMatriseA(A,n,h);

		//beregner de analytiske verdiene:
		double d = 2.0 / (h * h);
		double a = -1.0 / (h * h);
		for (int i = 0; i < n; i++) {
			eig_values(i) = d + 2 * a * cos(((i + 1) * acos(-1)) / (n + 1.0));
		}
		cout << "Analystisk verdi: " << endl;;
		cout << eig_values << endl;
	}
	if (z == 1) { //med potensial med 1 e
		cout << "Beregner egenverdiene med potensial med 1e:" << endl;
		fylleMatriseA_potensial(A, n, h, rho_0, potensial_1e);
	}
	if (z == 2) { //med potensial med 2 e
		cout << "Beregner egenverdiene med potensial med 2e:" << endl;
		fylleMatriseA_potensial(A, n, h, rho_0, potensial_2e);
	}

	double max_offdiag = maxoffdiag(A, &k, &l, n);
	while (fabs(max_offdiag) > epsilon && (double)iterations < max_number_iterations) {
		max_offdiag = maxoffdiag(A, &k, &l, n);
		rotate(A, R, k, l, n);
		iterations++;
	}
	cout << "Antall av interasjoner: " << iterations << "\n";
	cout << "Nummerisk beregnet egenverdiene: " << endl;
	for (int i = 0; i < n; i++) {
		num_eig(i) =  A(i, i);
	}
	//sorterer elementene i stigende rekkefolge
	for (int i = 0; i < n; i++){
		for (int j = i + 1; j < n; j++){
			if (num_eig[j] < num_eig[i])
			{
				double temp = num_eig[i];
				num_eig[i] = num_eig[j];
				num_eig[j] = temp;
			}
		}
		cout << num_eig[i] << endl;
	}
	//
	relativError(num_eig, a_verdi);

	string fileout = "fil_400_";
	string argument = to_string(rho_N);
	fileout.append(argument);
	ofile.open(fileout);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	for (int z = 0; z < 4; z++) {
		ofile << setw(15) << setprecision(8) << relativError(num_eig, a_verdi)[z] << endl;
	}
	ofile.close();


	/*
	//bruker Arma til å finne ekstakte egenverdier:
	eig_sym(arma_eigValue, A);
	cout <<"Armadillo sin verdi: " << endl;
	cout <<arma_eigValue << endl;*/
	return 0;
}

void fylleMatriseA(mat& A,int n, double h){
	//putter inn A elementer
	for (int i = 0; i < n; i++) {
		A(i, i) = 2.0 / (h * h);
		if (i != n - 1) {
			A(i + 1, i) = -1.0 / (h * h);
			A(i, i + 1) = -1.0 / (h * h);
		}
	}
}
void fylleMatriseA_potensial(mat& A, int n, double h, int rho_0, vec& p){
	//putter inn A elementer med potensial
	for (int i = 0; i < n; i++) {
		A(i, i) = 2.0 / (h * h) + p(i); //adder til potensial
		if (i != n - 1) {
			A(i + 1, i) = -1.0 / (h * h);
			A(i, i + 1) = -1.0 / (h * h);
		}
	}
}
double maxoffdiag(mat& A, int* k, int* l,  int n){
	double max = 0.0;
	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			if (fabs(A(i,j)) > max) {
				max = fabs(A(i, j));
				*l = i;
				*k = j;
			}
		}
	}
	return max;
}
void rotate(mat& A, mat& R, int k, int l, int n)
{
	double s, c;
	if (A(k, l) != 0.0) {
		double t, tau;
		tau = (A(l, l) - A(k, k)) / (2 * A(k, l));
		if (tau > 0) {
			t = 1.0 / (tau + sqrt(1.0 + tau * tau));
		}
		else {
			t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
		}
		c = 1 / sqrt(1 + t * t);
		s = c * t;
	}
	else {
		c = 1.0;
		s = 0.0;
	}
	double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
	a_kk = A(k, k);
	a_ll = A(l, l);
	// forandrer elementene i matrisene ved bruk av ligninger fra Lecture notes.
	A(k, k) = c * c * a_kk - 2.0 * c * s * A(k, l) + s * s * a_ll;
	A(l, l) = s * s * a_kk + 2.0 * c * s * A(k, l) + c * c * a_ll;
	A(k, l) = 0.0;
	A(l, k) = 0.0;
	//legger til resten av elementene
	for (int i = 0; i < n; i++) {
		if (i != k && i != l) {
			a_ik = A(i, k);
			a_il = A(i, l);
			A(i, k) = c * a_ik - s * a_il;
			A(k, i) = A(i, k);
			A(i, l) = c * a_il + s * a_ik;
			A(l, i) = A(i, l);
		}
		// Beregner egenvektorene:
		r_ik = R(i, k);
		r_il = R(i, l);
		R(i, k) = c * r_ik - s * r_il;
		R(i, l) = c * r_il + s * r_ik;
	}
}
double* relativError(vec& nummerical, int* exact){
	double *r = new double[4];
	for (int i = 0; i < 4; i++) {
		r[i] = abs((nummerical[i] - exact[i]) / (exact[i])); //relativ feil.
	}
	return r;
}
