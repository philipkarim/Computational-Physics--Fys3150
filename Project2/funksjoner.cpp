#include "prosjekt2funksjoner.h"

#include <armadillo>
#include <cmath>
using namespace std;
using namespace arma;




/*
Jacobi's method for finding eigenvalues
eigenvectors of the symetric matrix A.
The eigenvalues of A will be on the diagonal
of A, with eigenvalue i being A[i][i].
The j-th component of the i-th eigenvector
is stored in R[i][j].
A: input matrix (n x n)
R: empty matrix for eigenvectors (n x n)
n: dimention of matrices



//Tatt utgangspunkt i jacobialgoritme koden fra lectures2015

//maxxoffdiag funksjonen må endres på (kopiert)


*/



//finner max greiene
double maxoffdiag ( double ** A, int * k, int * l, int n ){
    double max = 0.0;
    for ( int i = 0; i < n; i++ ) {
        for ( int j = i + 1; j < n; j++ ) {
            if ( fabs(A[i][j]) > max ) {
    max = fabs(A[i][j]);
    *l = i;
    *k = j;
}}}
    return max;}


// Function to find the values of cos and sin
void rotate ( mat& A, mat& C, int k, int l, int n ){
    double s, c;
    if ( A[k][l] != 0.0 ) {
        double t, tau;
        tau = (A[l][l] - A[k][k])/(2*A[k][l]);
        if ( tau > 0 ) {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));}
        else {
            t = -1.0/( -tau + sqrt(1.0 + tau*tau));}

        c = 1/sqrt(1+t*t);
        s = c*t;
        }
        else {
        c = 1.0;
        s = 0.0;}

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];
    // changing the matrix elements with indices k and l
    A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
    A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
    A[k][l] = 0.0; // hard-coding of the zeros
    A[l][k] = 0.0;
// and then we change the remaining elements
    for ( int i = 0; i < n; i++ ) {
        if ( i != k && i != l ) {
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = c*a_ik - s*a_il;
            A[k][i] = A[i][k];
            A[i][l] = c*a_il + s*a_ik;
            A[l][i] = A[i][l];}
// Finally, we compute the new eigenvectors
        r_ik = C[i][k];
        r_il = C[i][l];
        C[i][k] = c*r_ik - s*r_il;
        C[i][l] = c*r_il + s*r_ik;
    }
    return;
}
//selve metoden
void jacobi_method (mat& A, mat& C, int n){
    // Setting up the eigenvector matrix
    for ( int i = 0; i < n; i++ ) {
        for ( int j = 0; j < n; j++ ) {
            if ( i == j ) {
                C[i][j] = 1.0;}
            else {
                C[i][j] = 0.0;}}}

    int k, l;
    double epsilon = 1.0e-8;
    double max_number_iterations = (double) n * (double) n * (double) n;
    int iterations = 0;
    double max_offdiag = maxoffdiag ( A, &k, &l, n );


    while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ){
        max:offdiag = maxoffdiag ( A, &k, &l, n );
        rotate ( A, C, k, l, n );
        iterations++;}

    cout << "Number of iterations: " << iterations << "\n";

    return;}
