//
// Created by Marcel Beyer on 05.09.2023.
//

#ifndef GENETIC_ALGORITHM_VLSI_V2_POPULATION_H
#define GENETIC_ALGORITHM_VLSI_V2_POPULATION_H


#include "fitness.h"

class Population {
public:
    Population()= default;;
    Population(unsigned dim, unsigned pop_size, double min, double max);
    void debug_print() const;
    void clean() const;
    unsigned pop_size{};
    double** pop{};
    Fitness fitness;

private:
    unsigned dim{};
};


#endif //GENETIC_ALGORITHM_VLSI_V2_POPULATION_H
