#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H
#include "planet.h"
#include <vector>

using std::vector;

class solarsystem
{
public:

    int number_planets;

    solarsystem();
    vector<planet> all_planets;
    void add(planet n);
    void print_position(vector<planet> vec);
    void print_position(vector<planet> vec, int n);

    void generate_matrixs(vector<planet> vec);
    void synctroniz(vector<planet> vec, double **mat);
};


#endif // SOLARSYSTEM_H
