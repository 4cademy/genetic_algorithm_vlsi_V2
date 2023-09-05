//
// Created by Marcel Beyer on 05.09.2023.
//

#include <random>
#include "population.h"

extern unsigned dim;
extern std::mt19937 gen;
extern std::uniform_real_distribution<double> uniformRealDistribution;
double** pop;
unsigned pop_size;


Population::Population(unsigned pop_size, double min, double max) {
    this->pop_size = pop_size;
    pop = new double*[pop_size];
    for (unsigned i = 0; i < pop_size; i++) {
        pop[i] = new double[dim];
    }
    for (unsigned i=0; i<pop_size; i++) {
        for (unsigned j=0; j<dim; j++) {
            pop[i][j] = uniformRealDistribution(gen) * (max-min) + min;
        }
    }
}

void Population::debug_print() const {
    for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < pop_size; j++) {
            printf("%f \t", pop[j][i]);
        }
        printf("\n");
    }
}

void Population::clean() const {
    for (unsigned i = 0; i < pop_size; i++) {
        delete[] pop[i];
    }
    delete[] pop;
}
