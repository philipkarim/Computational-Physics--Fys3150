//JACOBI OPPGAVE 2B
#include <iostream>
#include <cmath>
#include <armadillo>
#include <vector>
#include <array>

using namespace arma;
using namespace std;

double maxoffdiag(mat& A, int* k, int* l, int n);
void rotate(mat& A, mat& R, int k, int l, int n);

int main(){
    int n = 5;
    int N = n;
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
    cout << max_offdiag << endl;
    cout << R << endl;
    cout << A << endl;
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

    //Tester: Vi tester at funksjonen gir max og at funksjonen gir riktig egenverdi



    mat Atest = zeros<mat>(N, N);
    mat Ctest = eye<mat>(N, N);

    int * ar = new int[N];
    ar[0] = 1;
    ar[1] = 2;
    ar[2] = 3;
    ar[3] = 4;
    ar[4] = 5;
/*
    for(int i = 0; i < N; i++ ){
        for (int j = 0; j < N; j++){
            Atest(i,j) = ar[i];
        }

    }
    std::vector<int> dest(ar, ar + N);

        for (int i: dest) {
            std::cout << i << " ";
        }
    */

    for(int i = 0; i < N; i++ ){
        for (int j = 0; j < N; j++){
            if (abs(i-j)==1){Atest(i,j)=ar[1];}
            else if (abs(i-j)==2){Atest(i,j)=ar[2];}
            else if (abs(i-j)==3){Atest(i,j)=ar[3];}
            else if (abs(i-j)==4){Atest(i,j)=ar[4];}
            else {Atest(i,j)==ar[0];}
               }
            }
    cout << Atest << endl;
    //vec A = randu<vec>(5);
    //mat Atest = toeplitz();
    //setter max verdien vår til 9
    //Atest(3,2)=90;
    cout <<Atest<< endl;

    //maxoffdiag(Atest,&k,&l,N);

    //int max=9;
    cout << "Max er:  "<< maxoffdiag(Atest,&k,&l,N) << endl;
    if ((abs((maxoffdiag(Atest,&k,&l,N))-5.0)<epsilon)){
        cout<< "Test bestatt1"<<endl;}
    else {
        cout<<"Test ikke bestatt1"<<endl;
    }

//Riktige egenverdier test

    double * ev_jacobi_test = new double[N];
    double * ev_jacobi_test2 = new double[N];

    double * ev_fasit_test=new double[N];

    rotate(Atest,Ctest,k,l,N);

    for(int i = 0; i < N; i++ ){
            Ctest(i,i) = ev_jacobi_test[i];
        }

    double temp;


    //Loop som sorterer
    for(int i=0; i<N; i++)
        {
        for(int j=i+1; j<N; j++)
            {
                if(ev_jacobi_test[j] < ev_jacobi_test[i])
                {
                    temp = ev_jacobi_test[i];
                    ev_jacobi_test[i] = ev_jacobi_test[j];
                    ev_jacobi_test[j] = temp;
                }}
            }


    ev_fasit_test[0]=-0.94097;
    ev_fasit_test[1]=0;
    ev_fasit_test[2]=0;
    ev_fasit_test[3]=0;
    ev_fasit_test[4]=15.940971508067074;

    if (abs(ev_fasit_test[0]-ev_jacobi_test[0])<epsilon &&
        abs(ev_fasit_test[1]-ev_jacobi_test[1])<epsilon &&
        abs(ev_fasit_test[2]-ev_jacobi_test[2])<epsilon &&
        abs(ev_fasit_test[3]-ev_jacobi_test[3])<epsilon &&
        abs(ev_fasit_test[4]-ev_jacobi_test[4])<epsilon){
        cout<<"Test bestatt2"<<endl;}
    else {
        cout<<"Test ikke bestatt2"<<endl;
    }


    return 0;
}

double maxoffdiag(mat& A, int* k, int* l, int n){
    double max = 0.0;
    for (int i = 0; i < n; i++){
        for (int j = i+1; j < n; j++){
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
