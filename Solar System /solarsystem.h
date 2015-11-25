#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H
#include <armadillo>
#include "planet.h"
#include <vector>

using std::vector;

class solarsystem
{
public:

    int number_planets=0;

    solarsystem();
    vector<planet> all_planets;
    void add(planet n);
    void print_position(vector<planet> vec);
    void print_position(vector<planet> vec, int n);
    void synctroniz(vector<planet> vec, arma::mat &ma);
    void insert_data(vector<planet> vec, arma::mat &ma);

};


#endif // SOLARSYSTEM_H
