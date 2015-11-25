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

    double**  y_i     = new double*[7];
    double** y_i_temp = new double*[7];
    double** k1= new double*[7];
    double** k2 = new double*[7];
    double** k3 = new double*[7];
    double** k4 = new double*[7];

    solarsystem();
    vector<planet> all_planets;
    void add(planet n);
    void print_position(vector<planet> vec);
    void print_position(vector<planet> vec, int n);
    void synctroniz(vector<planet> vec, double **ma);
    void insert_data(vector<planet> vec, double **ma);

};


#endif // SOLARSYSTEM_H
