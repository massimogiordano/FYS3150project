#include <iostream>
#include <solarsystem.h>
#include <planet.h>
#include <cmath>
#include <armadillo>

using namespace arma;
using namespace std;

void sum_matrix(double **result, int coeff_one, double **first,int coeff_two, double **second, int n){
    for(int i=0; i<7; i++){
         for(int j=0; j<n; j++){
            result[j][i] = coeff_one*first[j][i] + coeff_two*second[j][i];
         }
    }
}
void printmat(double **ma, int n){
    cout << endl;
    for(int i=0; i<7; i++){

        for(int k=0; k<n;k++){
            cout <<  ma[k][i]<<" " ;
        } cout << endl;}
}
double force(double x, double y, double z, double Mothers){
    double G=6.67e-11;
    double force=0;
    double distance=0;

    distance = x*x + y*y + z*z;

    force = G*Mothers/pow(distance, 1.5);

    return force;
}
void derivate(double **dat, double **de, int n){

    double accelleration_x=0,accelleration_y=0,accelleration_z=0, mod_force;
//cout << endl;
    for(int i=0; i<n; i++){

        accelleration_x=0,accelleration_y=0,accelleration_z=0;
        for(int j=0; j<n; j++){
            if(i!=j){

mod_force = force(dat[j][0]-dat[i][0],dat[j][1]-dat[i][1] ,dat[j][2]-dat[i][2],  dat[j][6]);


accelleration_x += mod_force*(dat[j][0]-dat[i][0]);
accelleration_y += mod_force*(dat[j][1]-dat[i][1]);
accelleration_z += mod_force*(dat[j][2]-dat[i][2]);
}
}
        de[i][3] = accelleration_x; //velx
        de[i][4] = accelleration_y; //vely
 // **       cout << accelleration_y << "<- acc y"<< endl;
 // **       cout << accelleration_x << " " << de[i][3] <<"<- acc x"<< endl;
        de[i][5] = accelleration_z; //velz

    }

//cout<< "accelerazione: " << accelleration_x <<" "<< accelleration_y <<" "<< accelleration_z<< endl;

    for(int i=0; i<n; i++){
        de[i][0] = dat[i][3]; //velx
        de[i][1] = dat[i][4]; //velyjkjhjk
        de[i][2] = dat[i][5]; //vELz

    }
}


int main()
{

    solarsystem mysystem;

    //add a planet    Mass, x,y,z, vx,vy,vz;
    planet sun(   1.98e30, 0, 0,     0,0,0,   0 );
    planet earth( 5.24e24, 1.5e11, 0 ,0,0,30000,0);

    planet marth( 6.41e23, 2.2e11 ,0,0,0,24000,0);
    planet giove(1.89e27,7.78e11, 0,0,0,13000,0);
    planet saturno(5.7e26,1.4e12,0,0,0,9000,0);
    planet k( 6.41e23, 0, 2.2e1,0,0,24000,2);
    planet satellite(5.24e2, -1.5e11,0,0,0,30000,0);
    planet luna(7.3e2, 1.5e11, 3.84e5,0,-1000, 30000,0);

    mysystem.add(sun);
    mysystem.add(earth);
    //mysystem.add(marth);
    //mysystem.add(giove);
    //mysystem.add(saturno);
    //mysystem.add(satellite);
    //mysystem.add(luna);


    cout << sun.mass;
    mysystem.print_position(mysystem.all_planets,3);

    int elements = mysystem.number_planets;
    cout << "number of element" << elements<< endl;

    double**  y_i     = new double*[7];
    double** y_i_temp = new double*[7];
    double** k1= new double*[7];
    double** k2 = new double*[7];
    double** k3 = new double*[7];
    double** k4 = new double*[7];

    for(int i=0; i<7 ; i++){
        y_i[i] = new double[elements];
        y_i_temp[i] = new double[elements];
        k1[i] = new double[elements];
        k2[i] = new double[elements] ;
        k3[i] = new double[elements];
        k4[i] = new double[elements];
    }

    mysystem.insert_data(mysystem.all_planets, y_i);

cout <<"memorizzazione matrice:"<< endl;

printmat(y_i, elements);

double t=0, h=500, tmax=6e7;

cout << "succede qualcosa??"<< endl;



   while(t<tmax){

        derivate(y_i, k1, elements);



        sum_matrix(y_i_temp, 1, y_i, 0.5*h, k1, elements);
        derivate(y_i_temp, k2, elements);



        sum_matrix( y_i_temp, 1,  y_i, 0.5*h,  k2, elements);

        derivate( y_i_temp,  k3, elements);

        sum_matrix( y_i_temp, 1,  y_i, h,  k3, elements);

        derivate( y_i_temp,  k4, elements);
             for(int i=0; i<5; i++){
        for(int j=0; j<elements; j++){
               y_i[j][i] = y_i[j][i] + h*(k1[j][i] + 2*k2[j][i] + 2*k3[j][i] + k4[j][i])/6;
             }
        }


for(int j=elements-1; j<elements; j++){
cout << y_i[j][0] << "               " << y_i[j][1] << "            " ;
cout << endl;
}

    mysystem.synctroniz(mysystem.all_planets, y_i);

  t+=h;
   }



}








