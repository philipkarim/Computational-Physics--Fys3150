//JACOBI OPPGAVE 2A
#include <iostream>
#include <cmath>
#include <armadillo>
using namespace arma;
using namespace std;

double maxoffdiag(mat& A, int* k, int* l, int n);
void rotate(mat& A, mat& R, int k, int l, int n);

int main(){
	int n = 5;
	int k, l;
	int rho_0 = 0; int rho_N = 1;

	double epsilon = 1e-10;
	double max_number_iterations = n * n * n;
	int iterations = 0;

	//defining the matrix A and putting in elements:
	mat R = zeros<mat>(n, n);
	mat A = zeros<mat>(n, n);
	vec eig_values = zeros<vec>(n);
	vec eigenvectors = zeros<vec>(n);
	vec arma_eigValue = zeros<vec>(n);


	double h = (rho_N - rho_0) / ((float)n);
	//beregner diag og ikke_d_value elementer
	double d_value = 2.0 / (h * h);
	double ikke_d_value = -1.0 / (h * h);

	//putter inn A elementer
	for (int i = 0; i < n; i++) {
		A(i, i) = d_value;
		if (i != n - 1) {
			A(i + 1, i) = ikke_d_value;
			A(i, i + 1) = ikke_d_value;
		}
	}

	//putter inn R elementer
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				R(i,j) = 1.0;
			}
			else {
				R(i,j) = 0.0;
			}
		}
	}

	double max_offdiag = maxoffdiag(A, &k, &l, n);
	while (fabs(max_offdiag) > epsilon && (double)iterations < max_number_iterations) {
		max_offdiag = maxoffdiag(A, &k, &l, n);
		rotate(A, R, k, l, n);
		iterations++;
	}
	cout << "Number of iterations: " << iterations << "\n";
	cout << R << endl;
	cout << A << endl;

	//beregner de analytiske verdiene:
	double d = 2.0 / (h * h);
	double a = -1.0 / (h * h);
	for (int i = 0; i < n; i++) {
		eig_values(i) = d + 2 * a * cos(((i + 1) * acos(-1)) / (n + 1.0));
	}
	cout << eig_values << endl;
	//bruker Arma til å finne ekstakte egenverdier:
	eig_sym(arma_eigValue, A);
	cout << arma_eigValue << endl;

	return 0;
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
