
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
//#include <lib.h>

#define N 4



using namespace  std;
// output file as global variable
ofstream ofile;
// function declarations
void derivatives(double *, double *);
void initialise ( double&, double&, int&);
void output( double, double *, double);
void runge_kutta_4(double *, double *, int, double, double, double *, void(*)(double *, double *));

int main(int argc, char* argv[])
{
//  declarations of variables
  double *y, *dydt, *yout, t, h, tmax, E0;
  double initial_x, initial_v;
  int i, number_of_steps, n;
  char *outfilename;
  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] <<
      " read also output file on same line" << endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }
  ofile.open(outfilename);
  //  this is the number of differential equations
  n = 2;
  //  allocate space in memory for the arrays containing the derivatives
  dydt = new double[n];
  y = new double[n];
  yout = new double[n];
  // read in the initial position, velocity and number of steps
  initialise (initial_x, initial_v, number_of_steps);
  //  setting initial values, step size and max time tmax
  h = 4.*acos(-1.)/( (double) number_of_steps);   // the step size
  tmax = h*number_of_steps;               // the final time
  y[0] = initial_x;                       // initial position
  y[1] = initial_v;                       // initial velocity
  t=0.;                                   // initial time
  E0 = 0.5*y[0]*y[0]+0.5*y[1]*y[1];       // the initial total energy
  // now we start solving the differential equations using the RK4 method
  while (t <= tmax){
    derivatives(y, dydt);   // initial derivatives
    runge_kutta_4(y, dydt, n, t, h, yout, derivatives);
    for (i = 0; i < n; i++) {
      y[i] = yout[i];
    }
    t += h;
    output(t, y, E0);   // write to file
  }
  delete [] y; delete [] dydt; delete [] yout;
  ofile.close();  // close output file
  return 0;
}   //  End of main function

//     Read in from screen the number of steps,



//     initial position and initial speed
void initialise (double& initial_x, double& initial_v, int& number_of_steps)
{
 cout << "Initial position = ";
 cin >> initial_x;
 cout << "Initial speed = ";
 cin >> initial_v;
 cout << "Number of steps = ";
 cin >> number_of_steps;
}  // end of function initialise




//   this function sets up the derivatives for this special case
void derivatives(double *y, double *dydt)
{

   for(int i=0; i<N ; i+=2){
  dydt[0]=y[1];    // derivative of x
  dydt[1]=; // derivative of v
   }
} // end of function derivatives
//perchè passiamo il tempo se non lo usiamo??





//    function to write out the final results
void output(double t, double *y, double E0)
{

  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << t;
  ofile << setw(15) << setprecision(8) << y[0];
  cout << y[0]<< endl;
  ofile << setw(15) << setprecision(8) << y[1];
  ofile << setw(15) << setprecision(8) << cos(t);
  ofile << setw(15) << setprecision(8) <<
    0.5*y[0]*y[0]+0.5*y[1]*y[1]-E0 << endl;
}  // end of function output

/*   This function upgrades a function y (input as a pointer)
     and returns the result yout, also as a pointer. Note that
     these variables are declared as arrays.  It also receives as
     input the starting value for the derivatives in the pointer
     dydx. It receives also the variable n which represents the
     number of differential equations, the step size h and
     the initial value of x. It receives also the name of the
     function *derivs where the given derivative is computed
 */
void runge_kutta_4(double *y, double *dydx, int n, double x, double h,
                   double *yout, void (*derivs)( double *, double *))
{
  int i;
  double      xh,h_2,h6;
  double *dym, *dyt, *yt;
  //   allocate space for local vectors

  dym =  new double [n];
  dyt =  new double [n];
  yt  =  new double [n];

  h_2 = h*0.5;
  h6  = h/6.;
  xh  = x+h_2;

  for (i = 0; i < n; i++) {
    yt[i] = y[i]+h_2*dydx[i];
  }

  (*derivs)(yt,dyt);     // computation of k2, eq. 3.60

  for (i = 0; i < n; i++) {
    yt[i] = y[i]+h_2*dyt[i];
  }

  (*derivs)(yt,dym); //  computation of k3, eq. 3.61
  for (i=0; i < n; i++) {
     yt[i] = y[i]+h*dym[i];
    dym[i] += dyt[i];
  }
  (*derivs)(yt,dyt);    // computation of k4, eq. 3.62
  //      now we upgrade y in the array yout

  for (i = 0; i < n; i++){
    yout[i] = y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
  }

  delete []dym;
  delete [] dyt;
  delete [] yt;
}

//  end of function Runge-kutta 4
