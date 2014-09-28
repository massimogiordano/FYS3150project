#include <iostream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <time.h>

using namespace arma;
using namespace std;

#define NUMBER_OF_ELECTRON 1   //set 1 or 2
#define COULOMB_INTERACTION 1  //set 0 for NOT consender coulomb interaction
#define TOLLERANCE 10e-8
#define OMEGA 1
#define STEPS 100

double potential(double);
void  print_eigen_vector(mat &, int);
void biggest_element_of_matrix(mat &M, int *a, int *b);
double jacobi_algorithm(mat &M, mat &EIGEN_VECTORS);
void print_file(vec &obtained, vec &compared, double time);

int main()
{
    int n;
    double h, one_h2, p_i, temp=0, p_max=5;

    n=STEPS;

    if(NUMBER_OF_ELECTRON == 2) p_max=5/sqrt(OMEGA);

    h = p_max/double(n);
    one_h2 = 1/(h*h);

    cout << "Down and up diagonale values: " << one_h2 << endl;

    n--;                                               //matrix has dimension n - 1

    mat M(n,n);
    mat EIGENVECTOR(n,n);
    M.zeros();
    EIGENVECTOR.zeros();

    p_i = h;                                           //P_i is the distance from 0.
    M(0,0) = 2*one_h2 + p_i*p_i;
    EIGENVECTOR(0,0) = 1.;

    for(int i=1; i<n; i++){

        p_i = double(i+1)*h;

        M(i,i)    = 2.*one_h2 + potential(p_i);
        temp     +=   M(i,i);
        M(i-1, i) =  -one_h2;
        M(i, i-1) =  -one_h2;
        EIGENVECTOR(i,i) =1.;
    }

    cout << "Mediam value of diagonal elements: " <<  temp/n << endl;

    double time = jacobi_algorithm(M, EIGENVECTOR);

    cout << "Time need to diagolizzate the matrix: "<< time << "s" << endl;

    vec eigenvalue_jacobi(n);

    for(int i=0; i<n ; i++){             //it take the element on the diagonal (eigenvalue) and
       eigenvalue_jacobi(i) = M(i,i) ;   //it put it in a vector
    }

    for(int i=0; i<n ; i++){                     //This reorders the vector in ascending order.
        double temp =  eigenvalue_jacobi(i);
        double lowest = temp;

        int k=i;
        for(int j=i+1; j<n; j++){
            if( eigenvalue_jacobi(j) < lowest ){
            lowest = eigenvalue_jacobi(j);
            k = j;
            }
        }

        eigenvalue_jacobi(i) = lowest;
        eigenvalue_jacobi(k) = temp;

        for(int e=0; e<n ; e++){	       //This reorders the eigen vectors in the same order
            temp = EIGENVECTOR(e,i);
            EIGENVECTOR(e,i) = EIGENVECTOR(e,k);
            EIGENVECTOR(e,k) = temp;
        }
    }

    vec eigenvalue_armadillo(n);
    eigenvalue_armadillo = eig_sym(M);

    print_file(eigenvalue_jacobi, eigenvalue_armadillo, time);
    print_eigen_vector(EIGENVECTOR, 3);
}


double jacobi_algorithm(mat &M, mat &EIGEN_VECTORS){

    int n = M.n_rows;
    int iteration=0, break_loop=0;
    static char bar[] = "                                    "
                       "                                   ►";
    clock_t start= clock();

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
return double(clock() - start)/CLOCKS_PER_SEC;
}

double potential(double x){
         double potential=0;
         if(NUMBER_OF_ELECTRON == 2){
                potential = OMEGA*OMEGA*x*x;
                if (COULOMB_INTERACTION)
                    potential += 1/x;
         }else{
                potential = x*x;
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

//PRINT A FILE EIGEN VALUES
void print_file(vec &obtained, vec &compared, double time){
    int n = obtained.n_elem  ;
    char *filename = new char[1000];
    sprintf(filename, "Eigen_value_%d.txt", n+1);

        ofstream output (filename);
          if (output.is_open()){
            //TITLE
            output.precision(5);
            output << "Dimension matrix computed: "<< n << endl << "time for diagonalize: "<< time << "s"<< endl << endl;
            output <<std::scientific << "EIGEN VALUE" << "    Relative error " ;
            if(NUMBER_OF_ELECTRON==1)
                   output <<  "  Theoretical value   " << "   Relative error from teory: ";
            output << endl;

            for(int i=0; i<n ; i++){
                output <<std::scientific << obtained(i) << "    " << abs(compared(i)-obtained(i))*100/compared(i) << " %        ";
                if(NUMBER_OF_ELECTRON==1)
                       output <<  i*4 + 3 << "                    " << abs(3+ i*4-obtained(i))*100/(3+i*4) << " %";
                output << endl;
            }

            output.close();
          }
          else cout << "Unable to open file";
}
//PRINT A FILE EIGENVECTORS
void  print_eigen_vector(mat &EIGENVECTOR, int r){
    int n = EIGENVECTOR.n_rows;
    char *filename = new char[1000];
    sprintf(filename, "%d_Eigen_vectors_%d.txt",r, n+1);

        ofstream output (filename);
          if (output.is_open()){
            output.precision(5);
            output << "Dimension matrix computed: "<< n << endl  << endl;
            output <<std::scientific << "EIGEN VALUE"<< endl;

            for(int i=0; i<n ; i++){
                for(int y=0; y<r ; y++)
                output << EIGENVECTOR(i,y) << "     ";
                output << endl;
            }

            output.close();
          }
          else cout << "Unable to open file";
}
