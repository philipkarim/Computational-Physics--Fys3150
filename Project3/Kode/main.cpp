#include <iostream>
#include <cmath>

#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

using namespace std;


double uttrykk1(double x1, double y1, double z1, double x2, double y2, double z2);
double uttrykk2(double the1, double the2, double r_lag1, double r_lag2, double phi1, double phi2);
void gauleg(double, double, double *, double *, int);
double gammln(double);
void gaulag(double *, double *, int);

double pi=3.141592653, maxi=10.0, mini=1E-10; ;
int alp=2;

//Antall iterasjoner og integrasjonsgrenser
int n=10;
double x1=-5;
double x2=5;
//Definerer arrays og vektene til laguerre og legendre
double *xx = new double[n];
double *wx = new double[n];
double *rad = new double [n + 1];
double *wrad = new double [n + 1];
double *thet = new double [n];
double *wthet = new double [n];
double *xphi = new double [n];
double *wphi = new double [n];

//Korrekt svar
double exact = (5*pow(pi,2))/(pow(16.0,2));

int main(){
    //   loops for computations
    double int_gauss_leg;
    double int_gauss_lag;

    gauleg(x1,x2, xx, wx, n);
    gauleg(0,pi,thet,wthet,n);
    gauleg(0,2*pi,xphi,wphi,n);
    gaulag(rad,wrad,n);

    //Legendre løkke
    for (int i = 0;  i < n; i++){
        for (int j = 0; j < n; j++){
            for (int k = 0; k < n; k++){
                for (int l = 0; l < n; l++){
                    for (int m = 0; m < n; m++){
                        for (int o = 0; o < n; o++){
                            double weights = wx[i]*wx[j]*wx[k]*wx[l]*wx[m]*wx[o];
                            int_gauss_leg += weights*uttrykk1(xx[i],xx[j],xx[k],xx[l],xx[m],xx[o]);
    }}}}}}
    //Laguerre løkke
    for ( int i = 1;  i <= n; i++){
        for (int j = 1; j <= n; j++){
            for (int k = 0; k < n; k++){
                for (int l = 0; l < n; l++){
                    for (int m = 0; m < n; m++){
                        for (int o = 0; o < n; o++){
                            double weights = wrad[i]*wrad[j]*wthet[k]*wthet[l]*wphi[m]*wphi[o];
                            int_gauss_lag += weights*uttrykk2(rad[i],rad[j],thet[k],thet[l],xphi[m],xphi[o]);
    }}}}}}

    //multipliserer med konstanten fra endring i variabler
    int_gauss_lag *= pow((2*alp),-5);

          cout << "Legendre = " << int_gauss_leg <<"|| Korrekt losning: "<<exact<< "|| Feil:  " << abs(exact-int_gauss_leg)/exact << endl;
          cout << "Laguerre = " << int_gauss_lag <<"|| Korrekt losning: "<<exact<< "|| Feil:  " << abs(exact-int_gauss_lag)/exact << endl;

}


    //Uttrykk som skal integreres i kartesiske koordinater
    double uttrykk1 (double x1, double y1, double z1, double x2, double y2, double z2){
        double r1, r2,e1, e2, r12a, funk_a;

        //Definerer deler av uttrykket
        r1=sqrt(pow(x1,2)+pow(y1,2)+pow(z1,2));
        r2=sqrt(pow(x2,2)+pow(y2,2)+pow(z2,2));
        e1=exp(-2*alp*r1);
        e2=exp(-2*alp*r2);

        //Avstanden r12a holdes over 0 ved en if setning
        r12a=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
        funk_a=(e1*e2)/r12a;

        if(r12a < 1e-7){
            return 0;
        }
        else
            return funk_a;
    }

    //Uttrykk som skal integreres i kulekoordinater(Laguerre delen)
    double uttrykk2(double r_lag1, double r_lag2, double the1, double the2, double phi1, double phi2){
        double er1, er2, funk_b,funk_c, cosB, jaco, r12b;

        //Oppgitte verdier
        cosB=cos(the1)*cos(the2)+sin(the1)*sin(the2)*cos(phi1-phi2);
        r12b = pow(r_lag1,2)+pow(r_lag2,2)-2*r_lag1*r_lag2*cosB;

        //Uttrykket som skal integreres
        funk_b=(sin(the1)*sin(the2))/sqrt(r12b);

        //Avstanden r12b holdes over 0 ved en if setning
        if(r12b < 1e-6) {
            return 0;
        }
        else
            return funk_b;
    }


    //Legendre funksjon, utgangspunkt i eksempelkode fra github
    void gauleg(double a, double b, double x[], double w[], int n)
    {
       int         m,j,i;
       double      z1,z,xm,xl,pp,s,r,t;
       double      *x_low, *x_high, *w_low, *w_high;

       // roots are symmetric in the interval
       m  = (n + 1)/2;
       xm = 0.5 * (b + a);
       xl = 0.5 * (b - a);

       // pointer initialization
       x_low  = x;
       x_high = x + n - 1;
       w_low  = w;
       w_high = w + n - 1;

       // loops over desired roots
       for(i = 1; i <= m; i++) {
          z = cos(pi * (i - 0.25)/(n + 0.5));

          do {
             t =1.0;
             r =0.0;

         //Rekursjonsformel, s=> Lj+1(x), r=> Lj(x) and t=> Lj−1(x)
         for(j=1;j<=n;j++)
         {
            s=r, r=t;
            t = ((2.0*j-1.0)*z*r-(j-1.0)*s)/j;
         }


         /*t is now the desired Legrendre polynomial. Next compute ppp
           its derivative by standard relation involving also r,
           polynomial of one lower order.*/

         pp = n * (z * t - r)/(z * z - 1.0);
         z1 = z;
         // Newton's method
          z  = z1 - t/pp;
          } while(fabs(z - z1) > 1E-8);

          /* Scale the root to the desired interval and put in its symmetric
             counterpart. Compute the weight and its symmetric counterpart */

          *(x_low++)  = xm - xl * z;
          *(x_high--) = xm + xl * z;
          *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
          *(w_high--) = *(w_low++);
       }
    }

    //Gaus laguerre med funksjon "gammln" funksjon hentet fra eksempelkode i github
    void gaulag(double *x, double *w, int n)
    {
        int i,its,j;
        double ai;
        double t,r,s,pp,z,z1;

        for (i=1;i<=n;i++) {
            if (i == 1) {
                z=(1.0+alp)*(3.0+0.92*alp)/(1.0+2.4*n+1.8*alp);
            } else if (i == 2) {
                z += (15.0+6.25*alp)/(1.0+0.9*alp+2.5*n);
            } else {
                ai=i-2;
                z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alp/
                    (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alp);
            }
            for (its=1;its<=maxi;its++) {
                t=1.0;
                r=0.0;
                for (j=1;j<=n;j++) {
                    s=r;
                    r=t;
                    t=((2*j-1+alp-z)*r-(j-1+alp)*s)/j;
                }
                pp=(n*t-(n+alp)*r)/z;
                z1=z;
                z=z1-t/pp;
                if (fabs(z-z1) <= mini) break;
            }
            if (its > maxi) cout << "too many iterations in gaulag" << endl;
            x[i]=z;
            w[i] = -exp(gammln(alp+n)-gammln((double)n))/(pp*n*r);
        }
    }
    // end function gaulag
    double gammln( double xx)
    {
        double x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
            24.01409824083091,-1.231739572450155,
            0.1208650973866179e-2,-0.5395239384953e-5};
        int j;

        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
    }
