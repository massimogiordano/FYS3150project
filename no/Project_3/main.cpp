#include <iostream>
#include <solarsystem.h>
#include <planet.h>
#include <cmath>
using namespace std;

void sum_matrix(double **result, int coeff_one, double **first,int coeff_two, double **second, int n){
    for(int i=0; i<7; i++){
         for(int j=0; j<n; j++){

            result[j][i] = coeff_one*first[j][i] + coeff_two*second[j][i];

         }
    }
}


    /* force(direction in wich compute force, [........] ) */

    double force(double x, double y, double z, double Mown, double Mothers){
    double G=6.67e-11;
    double force=0;
    double distance=0;

    //cout << "Funziona?"<< x <<" "<< y << " "<< z << endl;

    distance = x*x + y*y + z*z;
                                //cout << distance << " "<< Mothers << endl;

    force = -G*Mothers/pow(distance, 1.5);
                               // cout << "froce"<<force << endl;
    return force;
}


void derivate(double **data, double **de, int n){

// function to compute the forces on each planet.
    double accelleration_x=0,accelleration_y=0,accelleration_z=0, mod_force;

    for(int i=0; i<n; i++){
        accelleration_x=0,accelleration_y=0,accelleration_z=0;
        for(int j=0; j<n; j++){
            if(i!=j){

mod_force = force(data[j][0]-data[i][0],data[j][2]-data[i][2] ,data[j][4]-data[i][4],  data[i][6],  data[j][6]);
cout << data[j][6] << " " << j << "<-massa pianeta" << endl;


          //  cout << "Modulo accelerazione:"<< i <<"-"<< j << " " << mod_force << endl;

            accelleration_x += mod_force*(data[j][0]-data[i][0]);
            accelleration_y += mod_force*(data[j][2]-data[i][2]);
            accelleration_z += mod_force*(data[j][4]-data[i][4]);
}
        }
           cout<< "accelerazione: " << accelleration_x <<" "<< accelleration_y << endl;
            de[i][1] = accelleration_x; //velx
            de[i][3] = accelleration_y; //vely
        //    de[i][5] = accelleration_z; //velz

    }


/*  This gives to the second matrix the value of derivate
 *  and register the value of mass too.
 */


    for(int i=0; i<n; i++){
        de[i][0] = data[i][1]; //velx
        de[i][2] = data[i][3]; //vely
      //  de[i][4] = data[i][5]; //velz
        de[i][6] = data[i][6]; //mass

    }
//calcolare forza su ogni pianeta su ogni direzione.
//aggiungere un elemento ad array per mettere anche massa pianeta

}

int main()
{
    solarsystem mysystem;

    //add a planet    Mass, x,y,z, vx,vy,vz;
    planet sun(   1.98e30, 0, 0,     0,0,0,   0 );
    planet earth( 5.24e24, 1.5e11, 0 ,0,0,30000,0);

    planet marth( 6.41e1, 0, 2.2e11,0,0,24000,0);
    planet k( 6.41e23, 0, 2.2e1,0,0,24000,2);

    mysystem.add(sun);
    mysystem.add(earth);
   // mysystem.add(marth);


    mysystem.print_position(mysystem.all_planets,3);

   int elements = mysystem.number_planets;

//__________________________iniziare rudda kutta
//cout<< elements;
   double t=0, h=100, tmax=3e7;

cout << "succede qualcosa??"<< endl;


   while(t<tmax){




        derivate(mysystem.y_i, mysystem.k1, elements);

/*        sum_matrix(mysystem.y_i_temp, 1, mysystem.y_i, 0.5*h, mysystem.k1, elements);

        derivate(mysystem.y_i_temp, mysystem.k2, elements);

        sum_matrix(mysystem.y_i_temp, 1, mysystem.y_i, 0.5*h, mysystem.k2, elements);

        derivate(mysystem.y_i_temp, mysystem.k3, elements);

        sum_matrix(mysystem.y_i_temp, 1, mysystem.y_i, h, mysystem.k3, elements);

        derivate(mysystem.y_i_temp, mysystem.k4, elements);
*/
        //sun all eximation

        for(int i=0; i<6; i++){
        for(int j=0; j<elements; j++){

 if(i==0) cout << j << " - " << "   " << mysystem.y_i[j][0] <<" - " << "   " << mysystem.y_i[j][6];

//mysystem.y_i[j][i] = mysystem.y_i[j][i];// + h*(mysystem.k1[j][i] + 2*mysystem.k2[j][i]+2*mysystem.k3[j][i] + mysystem.k4[j][i])/6 ;

if(i==0) cout << " - " << "   " << mysystem.k1[j][2] << "   " << mysystem.k1[j][6]<< endl;

             }
        }

mysystem.synctroniz(mysystem.all_planets, mysystem.y_i);

//mysystem.print_position(mysystem.all_planets, 2);
//funziona aggiorna dati posizione e velocità pianeti

        t+=h;
   }








for(int i=0; i<7; i++){

     for(int j=0; j<elements; j++){
      // cout << mysystem.y_i[j][i]<< " p ";
     }
 cout << endl;

}


cout << "ciao!!"<< endl;


    return 0;
}

