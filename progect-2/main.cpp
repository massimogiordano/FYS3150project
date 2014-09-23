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
#define TOLLERANCE 10e-8
#define OMEGA 10e10
#define STEPS 100


int main()
{
    int n, i,j;
    double h, one_h2, p_i, tollerance, temp=0, lowest;

    n=STEPS;

    h = P_MAX / double(n);
    cout << "h: " << h << endl;
    one_h2 = 1/(h*h);
    tollerance = TOLLERANCE;

    n--;    //matrix dimension is n-1

    mat M(n,n);
    mat E(n,n);
    M.zeros();
    E.zeros();
    p_i = h;

    cout << endl << "Down and up diagonale values: " << one_h2 << endl;


    M(0,0) = 2*one_h2 + p_i*p_i;
    E(0,0) = 1.;

    for(i=1; i<n; i++){

        p_i = double(i+1)*h;

        M(i,i)    = 2.*one_h2 + p_i*p_i;
        temp += M(i,i);
        E(i,i)    = 1.;
        M(i-1, i) = -one_h2;
        M(i, i-1) = -one_h2;

    }

    cout << "Mediam value of diagonal elements: " <<  temp /n << endl;

    vec eigen_arma(n);
    eigen_arma = eig_sym(M);

    double big_a, teta, c, s, c2,s2,cs, M_kk, M_kl, M_ll,tau, tan;
    int k, l, count=0, break_loop=0;                //elemento sul quale ruotare.

    //START LOOP

    while(count < 1000000 && break_loop==0) {		//impostare numero max sicurezza
        big_a=0;
        count++;

    for(i=0; i<n; i++){                            //cerco il valore più grande
    for(int j=i+1; j<n; j++){
        if( abs(M(i,j)) > big_a){
            big_a = abs(M(i,j));
            k=i;
            l=j;
        }
    }}

    cout << big_a << "  k-l    "<< k << " "<< l<< endl;

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
        cs = c*s;
        c2= c*c;
        s2 = s*s;

        M_kk = M(k,k); M_kl = M(k,l); M_ll = M(l,l);               //salvo elelementi della matrice che verrebero riscritti

        //sovascrivo matrice

        M(k,k) = c2*M_kk - 2*cs*M_kl + s2*M_ll;
        M(l,l) = s2*M_kk + 2*cs*M_kl + c2*M_ll;
        M(k,l) = M(l,k) = 0;

        for(i=0; i<n ; i++){
            if(i!=k && i!=l){
                temp= M(k,i);

                M(k,i) = M(i,k) = c*M(k,i) - s*M(l,i);
                M(l,i) = M(i,l) = s*temp + c*M(l,i);
            }

            temp = E(i,k);
            E(i,k) = c*temp - s*E(i,l);
            E(i,l) = c*E(i,l) + s*temp;

        }

        //




    }else{
        cout << endl << " Number of iteration:" <<  count    << endl;
        break_loop = 1; //stop while
    }

    }

    vec eigen_jac(n);
                                    //cout << E << endl; //______
    for(i=0; i<n ; i++){
       eigen_jac(i) = M(i,i) ;
    }


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

        for(int e=0; e<n ; e++){
            temp = E(e,i);
            E(e,i) = E(e,k);
            E(e,k) = temp;
        }



    }

 for(i=0; i<n ; i++) cout << eigen_jac(i) << "rel err: " << abs(eigen_arma(i)-eigen_jac(i))*100/eigen_arma(i) << "%  " <<  i*4 + 3 << "   rel err from teo: " << abs(3+ i*4-eigen_jac(i))*100/(3+i*4) << "%"<< endl;

 for(i=0; i<n ; i++) cout << i <<  "  "<< E(i,0) << " " << E(i,1) << " " << E(i,2)<< endl;;



// trovare autovalori con jacobi rotazione


}
