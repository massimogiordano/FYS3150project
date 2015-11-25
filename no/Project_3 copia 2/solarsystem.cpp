#include "solarsystem.h"
#include "planet.h"
#include <iostream>
using namespace std;
solarsystem::solarsystem()
{
}
void solarsystem::add(planet n){
         all_planets.push_back(n);  //questa funzione inserisce un elementp nella push_bac
 }
void solarsystem::print_position(vector<planet> vec){
    print_position(vec, 3);
}

void solarsystem::print_position(vector<planet> vec, int n){
    if(n>3 || n<=0) n=3;
    for(int i=0; i<vec.size(); i++){
        planet &questo = vec[i];
        std::cout << std::scientific;
        for(int j=0; j<n;j++){
            std::cout << questo.position[i] << " ";
        }
        std::cout << "    ";
       }
    std::cout << std::endl;
}

void solarsystem::generate_matrixs(vector<planet> vec){
    number_planets = vec.size();
   std::cout <<"!! "<< number_planets << "!!!!"<<std::endl;
   double** y_i = new double*[7];
   double** y_i_temp = new double*[7];
   double** k1 = new double*[7];
   double** k2 = new double*[7];
   double** k3 = new double*[7];
   double** k4 = new double*[7];
    for (int i = 0; i < 7; i++){
        y_i[i] = new double[number_planets];
        y_i_temp[i] = new double[number_planets];
        k1[i] = new double[number_planets];
        k2[i] = new double[number_planets];
        k3[i] = new double[number_planets];
        k4[i] = new double[number_planets];
    }

    for(int j=0; j<number_planets;j++){
        planet &questo = vec[j];
        y_i[j][6] = questo.mass;

        for (int i = 0; i < 3; i++){
            y_i[j][2*i] = questo.position[i];
            y_i[j][(2*i)+1] = questo.velocity[i];
            }
        }
cout << "verifica memorizzazione matrice:"<< endl;
    for (int i = 0; i < 7; i++){
        for (int j = 0; j < number_planets; j++){
       cout<<y_i[j][i] << "_";
}
   cout << endl;
    }

}

void solarsystem::synctroniz(vector<planet> vec, double **mat){
    int n = vec.size();

    for(int j=0; j<n;j++){
        planet &questo = vec[j];
    for (int i = 0; i < 3; ++i){
       questo.position[i] =  y_i[j][2*i];
       questo.velocity[i] = y_i[j][(2*i)+1];
    }
}
}
