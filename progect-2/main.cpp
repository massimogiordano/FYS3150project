#include <iostream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <time.h>


using namespace arma;
using namespace std;

#define P_MAX  5
#define TOLLERANCE 10e-10
#define PI 3.14
#define OMEGA 10e10
#define STEPS 50


int main()
{
    int n, i,j;
    double h, one_h2, p_i, w=OMEGA, ht, tollerance;

    n=STEPS;



    h = P_MAX / double(n);
    cout << "h: " << h << endl;
    one_h2 = 1/(h*h);
    //ht = 6.62606957*pow(10.,-34.);
    n--; //matrix dimension is n-1

    mat M(n,n);
    mat B(n,n);
    M.zeros();

    p_i = h;
    M(0,0) = 2*one_h2 + p_i*p_i;
    cout << endl << "Down and up diagonale values: " << one_h2 << endl;
    cout << endl << "Vpotential P_max: " << P_MAX*P_MAX << endl;

    for(i=1; i<n; i++){

        p_i = double(i+1.)*h;  //chiedere circa il +1 e se bisogna partire da zero o da uno
                                                                                            //cout << p_i << endl;

        M(i,i) = 2.*one_h2 + p_i*p_i;
        M(i-1, i) = -one_h2;
        M(i, i-1) = -one_h2;

    }

    B = M;

    //stampo matrice
    cout << "Diagonal values of matrix: " << endl;
    for(i=0; i<n; i++){
        cout << M(i,i) << "  " ;
    }

//cout << M << endl;
    tollerance = TOLLERANCE;

    //iterazione qui

    vec eigen_arma(n);
    eigen_arma = eig_sym(B);

    cout << eigen_arma << endl;
    cout << "||||| autovettori da armadillo |||||";

    double big_a, teta, c, s, c2,s2,cs, M_kk, M_kl, M_ll,tau, tan;
    int k, l, count=0, break_loop=0;           //elemento sul quale ruotare.

    //START LOOP

    while(count < 1000 && break_loop==0) {		//impostare numero max sicurezza
        big_a=0;
        count++;
                                                         // cout << count << endl;
    for(i=0; i<n; i++){ //cerco il valore più grande
    for(int j=i+1; j<n; j++){
                                                         //cout << "in" << i << " "<< j << endl;
                                                         //cout << M(i,j) << " " << big_a << endl;
                 if( abs(M(i,j)) > big_a){
                     big_a = abs(M(i,j));
                     k=i;
                     l=j;
                                                         //cout << big_a << "  k-l    "<< k << j<< endl;
                 }

    }}

    if(abs(M(k,l))>tollerance){
        //diagonalizzation

        tau = (M(l,l) - M(k,k))/(2*M(k,l));
        if(tau > 0){
            tan = 1./(tau + sqrt(1.0 + tau*tau));
        }else{
            tan = -1./(-tau + sqrt(1.0 + tau*tau));
        }
        c=1./sqrt(1+(tan*tan));
        s=tan/sqrt(1+tan*tan);

/*
        if(M(l,l)!=M(k,k)){
        teta=1;
        //   teta = 0.5* atan((2*M(k,l))/(M(l,l)-M(k,k)));
        }else{
        teta = acos(-1)/4;
        }					//FINE TROVARE TETA.
        //calcolo seni e coseni.
        c = cos(teta); s= sin(teta); cs= c*s; c2= c*c; s2= s*s;
*/
        M_kk = M(k,k); M_kl = M(k,l); M_ll = M(l,l);               //salvo elelementi della matrice che verrebero riscritti

        //sovascrivo matrice

        M(k,k) = c2*M_kk - 2*cs*M_kl + s2*M_ll;
        M(l,l) = s2*M_kk + 2*cs*M_kl + c2*M_ll;
        M(k,l) = M(l,k) = cos(2*teta)*M_kl + 0.5*sin(2*teta)*(M_kk - M_ll);//(c2 - s2)*M_kl + cs*(M_kk - M_ll);

        for(i=0; i<n ; i++){
            if(i!=k && i!=l){
                M(k,i) = M(i,k) = c*M(k,i) - s*M(l,i);
                M(l,i) = M(i,l) = s*M(k,i) + c*M(l,i);
            }
        }

                                                        //cout << "inizio" << k << " "<< l << endl;

    }else{
        cout << endl << "!!!!!Number of iteration:" <<  count    << endl;
        break_loop = 1; //stop while
    }

    }

    double less_error, error, lowest=10, temp;
    vec eigen_jac(n);

    //eigen_jac = diagvec(M,0);
  //  eigen_jac(n+1) = eigen_jac(0)
 for(i=0; i<n ; i++){
 eigen_jac(i) = M(i,i) ;
 }

cout << eigen_jac << endl << endl;
//riordino i vettori per essere comparati con gli autovettori

    for(i=0; i<n ; i++){
        lowest = temp = eigen_jac(i);
        k=i;
        for(j=i+1; j<n; j++){
            if( eigen_jac(j) < lowest ){
            lowest = eigen_jac(j);
            k = j;

            }
        }

        eigen_jac(i) = lowest;
        eigen_jac(k) = temp;


    }



cout << eigen_jac << endl;
    //riordino gli autovettori
/*
    for(i=0; i<n ; i++){
        less_error = 10000;
    for(int j=0; j<n ; j++){
        error = (abs(abs(eigen_arma(j)) - abs(M(i,i)) ) / abs(eigen_arma(j)));
        if (less_error > error){
                less_error = error;
                l=j;
    }
    }


        eigen_same_order(i) = eigen_arma(l);

    }


    //stampo errori relativ

    error = 0;
    for(int j=0; j<n ; j++){
        if( (abs(abs(eigen_same_order(j)) - abs(M(j,j)) ) / abs(eigen_same_order(j))) > error)
            error = (abs(abs(eigen_same_order(j)) - abs(M(j,j)) ) / abs(eigen_same_order(j)));
    }

    cout << endl << "Maximum relativ error of EIGENVALUE copared with armadillo one is:  " << error << "%." << endl;

    cout << endl << "ARMADILLO VALUES:     JACOBI'S VALUES" << endl;

    for(int j=0; j<n ; j++){
    cout << endl << eigen_same_order(j) <<  "            " << M(j,j) << endl;
    }

*/



// trovare autovalori con jacobi rotazione


}
