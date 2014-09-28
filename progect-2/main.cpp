#include <iostream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <time.h>
//#include <unistd.h>

using namespace arma;
using namespace std;

#define NUMBER_OF_ELECTRON 1   //set 1 or 2
#define COULOMB_INTERACTION 0  //set 0 for NOT consender coulomb interaction
#define TOLLERANCE 10e-8
#define OMEGA 1/(132.638)
#define STEPS 100

double potential(double x){
 double potential=0;
         if(NUMBER_OF_ELECTRON == 2){
                potential = OMEGA*OMEGA*x*x;
                if (COULOMB_INTERACTION)
                    potential += 1/x;
         }

return potential;
}

void biggest_element_of_matrix(mat &M, int *a, int *b){

    double biggest_element=0;
    int n = M.n_rows;
    for(int i=0; i<n; i++){

    for(int j=i+1; j<n; j++){
            if( abs(M(i,j)) > biggest_element){
            biggest_element = abs(M(i,j));
            *a=i;
            *b=j;

        }
    }}

}


void jacobi_algorithm(mat &M, mat &EIGEN_VECTORS){

    int n = M.n_rows;
    int iteration=0, break_loop=0;
    static char bar[] = "                                    "
                        "                                   ►";

    cout << endl <<"Digonalizzation of Matrix in progress:"<< endl;
    double first = log(abs(M(0,1))), normalizzation=-log(TOLLERANCE);   //value to normalizze the time need to the algorithm

    while(iteration < 5*n*n && break_loop==0) {	// "Iteration" provides to count the number of iteration.
                                                    // If for some reason the alghorithm doesn't converge the
        iteration++;                                // loop stop after n*n interations. () This avoid an
                                                    // infinitive loop.

        int k,l;
        biggest_element_of_matrix(M,&k,&l);         //this function returns the index largest element

        //this formula can estimate the percentage of completion of the diagonalizzation
        int loading= ((-log(abs(M(k,l))) + first)*100/(normalizzation+first));
        //this make a very nice loading bar.
        printf("\x1b[33m\033[40m" "  %d%c %s\r ",loading, 37, &bar[40-int(loading/2.2)]);  fflush(stdout);

        if( abs(M(k,l)) < TOLLERANCE)  break_loop=1;  //when the biggest element is less than tollerance the algoritm can stop

        double  tau = (M(l,l) - M(k,k))/(2*M(k,l));        //it found the angol that after rotation make the element A_kl=0
        double  tan;
        if(tau > 0){
           tan = 1./(tau + sqrt(1.0 + tau*tau));//-tau + sqrt(1.0 + tau*tau);
        }else{
           tan =  -1./(-tau + sqrt(1.0 + tau*tau));//-tau - sqrt(1.0 + tau*tau);//
        }
        double c= 1./sqrt(1+(tan*tan));
        double s= tan/sqrt(1+tan*tan);
        double cs = c*s, c2= c*c, s2 = s*s;                    //Its compute cos^2, sen^2, cos*sen.
        double M_kk = M(k,k), M_kl = M(k,l), M_ll = M(l,l);    //Its memorize the elements of the matrix are rewritten


        M(k,l) = M(l,k) = 0;               					// This is the result of the trasformation for similarity
        M(k,k) = c2*M_kk - 2*cs*M_kl + s2*M_ll;				//
        M(l,l) = s2*M_kk + 2*cs*M_kl + c2*M_ll;				//

        for(int i=0; i<n ; i++){
            if(i!=k && i!=l){
                double temp= M(k,i);
                M(k,i) = M(i,k) = c*M(k,i) - s*M(l,i);      //
                M(l,i) = M(i,l) = s*temp + c*M(l,i);
            }
            double temp = EIGEN_VECTORS(i,k);
            EIGEN_VECTORS(i,k) = c*temp - s*EIGEN_VECTORS(i,l);
            EIGEN_VECTORS(i,l) = c*EIGEN_VECTORS(i,l) + s*temp;
        }
    }
    printf("\x1b[31m\033[0m" "\nNumber of iteration: %d \n", iteration);
}

int main()
{
    int n;
    double h, one_h2, p_i, temp=0, p_max=5;

    n=STEPS;

    if(NUMBER_OF_ELECTRON == 2) p_max=5/sqrt(OMEGA);

    h = p_max/double(n);
    one_h2 = 1/(h*h);

    cout << "Down and up diagonale values: " << one_h2 << endl;

    n--;                                                    //matrix dimension is n-1

    mat M(n,n);
    mat E(n,n);
    M.zeros();
    E.zeros();

    p_i = h;                                           //P_i is the distance from 0.
    M(0,0) = 2*one_h2 + p_i*p_i;
    E(0,0) = 1.;

    for(int i=1; i<n; i++){

        p_i = double(i+1)*h;

        M(i,i)    = 2.*one_h2 + potential(p_i); OMEGA*OMEGA*p_i*p_i + 1/p_i;  //p_i*p_i; //
        temp += M(i,i);
        E(i,i)    = 1.;
        M(i-1, i) = -one_h2;
        M(i, i-1) = -one_h2;

    }

    cout << "Mediam value of diagonal elements: " <<  temp /n << endl;

    vec eigen_arma(n);
    eigen_arma = eig_sym(M);


    jacobi_algorithm(M, E);

    vec eigen_jac(n);

    for(int i=0; i<n ; i++){
       eigen_jac(i) = M(i,i) ;
    }


//riordino i vettori per essere comparati con gli autovettori

    for(int i=0; i<n ; i++){
        double temp =  eigen_jac(i),  lowest = temp;

        int k=i;
        for(int j=i+1; j<n; j++){
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

 for(int i=0; i<n ; i++) cout << eigen_jac(i) << "rel err: " << abs(eigen_arma(i)-eigen_jac(i))*100/eigen_arma(i) << "%  " <<  i*4 + 3 << "   rel err from teo: " << abs(3+ i*4-eigen_jac(i))*100/(3+i*4) << "%"<< endl;

 for(int i=0; i<n ; i++) cout << i <<  "  "<< E(i,0)*E(i,0) << " " << E(i,1)*E(i,1) << " " << E(i,2)*E(i,2)<< endl;;





}
