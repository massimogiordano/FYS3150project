#include <iostream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <stdio.h>
#include <time.h>


using namespace arma;
using namespace std;

#define FUNX  100*exp(-10*x)
#define ANALYTICSOLUT (1 - (1-exp(-10))*double(i)*h - exp(-10.*double(i)*h))


double func_h2(double, double);
double func(double);
void LU_decomposition(int, double, double* );



int main()
{

    unsigned int n=1, i, p, t, lu;
    double *f_i, *u_i, h2, h, *v_i, *e_i, *lu_sol, e_max=0, time, time_lu;
    char name[10];
    clock_t start, finish;

    cout << "Number of iteration of 10^n: ";
    cin >> p;

    for(t=1; t<=p; t++){

    n *=10;
    cout << "N= "<< n << endl;


    h = 1./(1.*n + 1.);
    h2 = h*h;

    f_i = new double[n+1];
    v_i = new double[n+1];
    u_i = new double[n+2];
    e_i = new double[n+1];
    lu_sol = new double[n];

    f_i[0] = 0.;

    //f_i
    for(i=1; i<=n; i++ ){
        f_i[i] = func_h2(i*h, h2);
        v_i[i] = ANALYTICSOLUT;          		      //analytic solution
    }

    //numerical solution
    start = clock();             		      //start of solving alghortm
    for(i=1; i<=n; i++ ){
        f_i[i] = f_i[i-1] + i*f_i[i];       	  		        // 3 flop
    }

    u_i[n+1] = 0.;        		//to be able to use the generical formula

    for(i=n; i>=1; i--){
        u_i[i] =( double(i)*u_i[i+1] + f_i[i] )/double(i+1);    // 4 flop
    }

    finish = clock();
    time = double(finish-start)/CLOCKS_PER_SEC;
    cout <<"    numerical  "<< time << endl;

    //relative error
    for(i=1; i<=n; i++){
        e_i[i] = log10(abs(u_i[i]/ANALYTICSOLUT -1));
        if(e_i[i]<e_max){
               e_max=e_i[i];        //finding the maximum realative error
        }
    }


    if(n<100000){   //If the matrix is too big is better don't compute it
        lu = 0;              //set value to control the printing of value
        start=clock();
        LU_decomposition(n,h, lu_sol);
        finish=clock();
        time_lu = double(finish-start)/CLOCKS_PER_SEC;
        cout<<"    LU         "<< time_lu <<endl;
    }else{
        cout<<"    LU         "<< "No computed." <<endl;
        lu = 1;              //set value to control the printing of values
    }                        //zero means computed, one means Not computed



    name[0] = '1';	          //try to find a way to print the file called
    for(i=1; i<t+1;i++){      //as the numbers n used to find the solution
        name[i]='0';
    }

    name[i++]='.';
    name[i++]='t';
    name[i++]='x';
    name[i++]='t';
    name[i++]='\0';


    char time_lu_char[32];
    snprintf(time_lu_char, sizeof(time_lu_char), "%g", time_lu);

    ofstream output (name);
      if (output.is_open())
      {
        output.precision(5);
        output << "Dimension matrix computed: 10^" <<  t << "x10^" << t << endl << endl;
        output << "Duration of Numerical calculation:  " << time << "\n                   LU calculation:  " << (lu ? "NA" : time_lu_char)<< "\n \nMaximun relative error: " <<  e_max << endl;
        output << "x_i            u_i            LU solution    Analytic sol   f_i            relative error (Log10 Ei)\n\n";

        for(i=1; i<=n; i++){
            double ih = double(i)*h;

            output << std::scientific << ih << "    " << u_i[i] << "    " << lu_sol[i-1] << "    " << v_i[i] << "    " << func(ih) << "    " << e_i[i] << std::scientific<< endl;
        }

        output.close();
      }
      else cout << "Unable to open file";

    }
    return 0;
}


double func(double x){

    double value;
    value = FUNX;
    return value;
}


double func_h2(double x, double h2){

    double value;
    value = h2*FUNX;
    return value;
}

void LU_decomposition(int n, double h, double *lu_sol){
    vec f(n);
    mat M(n,n);
    M.zeros();
    M(0,0)=2;
    f[0]= func_h2(h, h*h);

    for(int i=1;i<n;i++){
        f[i]= func_h2((i+1.)*h, h*h);
        M(i,i)=2;
        M(i-1,i)=-1;
        M(i,i-1)=-1;
    }
    vec x=solve(M,f);

    for(int i=0; i<n; i++){
    lu_sol[i] = x[i];
    }
}
