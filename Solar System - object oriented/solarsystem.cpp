#include "solarsystem.h"
#include "planet.h"
#include <iostream>
#include <armadillo>
using namespace arma;
using namespace std;
solarsystem::solarsystem()
{
}
void solarsystem::add(planet n){
     number_planets++;
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
        std::cout << questo.position[j] << " ";
        }
        std::cout << "    ";
       }
    std::cout << std::endl;
}

void solarsystem::insert_data(vector<planet> vec, mat &ma){
    for(int i=0; i<vec.size(); i++){
        planet &questo = vec[i];
        ma(i,6)=questo.mass;

        for(int k=0; k<3;k++){
            ma(i,k)=questo.position[k];
            ma(i,k+3)=questo.velocity[k];
        }

    }



}

void solarsystem::synctroniz(vector<planet> vec, mat &ma){
    int n = vec.size();

    for(int j=0; j<n;j++){
        planet &questo = vec[j];
    for (int i = 0; i < 3; ++i){
       questo.position[i] =  ma(j,i);
       questo.velocity[i] = ma(j,i+3);
}
}

}
