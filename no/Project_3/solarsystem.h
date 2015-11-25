#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H
#include "planet.h"
#include <vector>

using std::vector;

class solarsystem
{
public:
    double** y_i = new double*[7];
    double** y_i_temp = new double*[7];
    double** k1 = new double*[7];
    double** k2 = new double*[7];
    double** k3 = new double*[7];
    double** k4 = new double*[7];
    int number_planets;

    solarsystem();
    vector<planet> all_planets;
    void add(planet n);
    void print_position(vector<planet> vec);
    void print_position(vector<planet> vec, int n);

    void derivate(vector<planet> vec, double **mat );
    void generate_matrixs(vector<planet> vec);
    void synctroniz(vector<planet> vec, double **mat);
};


#endif // SOLARSYSTEM_H
