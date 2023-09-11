//
// Created by Marcel Beyer on 05.09.2023.
//

#ifndef GENETIC_ALGORITHM_VLSI_V2_GA_OLD_H
#define GENETIC_ALGORITHM_VLSI_V2_GA_OLD_H


#include "population.h"

class GaOld {
public:
    GaOld(Population &population, unsigned dim);

    void evolve(int generations);
    void clean() const;

private:
    unsigned dim;
    unsigned pop_size;
    Population population;

    void selection_roulette() const;
    void crossover_uniform() const;
};


#endif //GENETIC_ALGORITHM_VLSI_V2_GA_OLD_H
