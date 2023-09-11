//
// Created by Marcel Beyer on 05.09.2023.
//

#include <random>
#include "ga_old.h"
#include "omp.h"

extern std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
extern std::uniform_real_distribution<double> uniformRealDistribution;
double** mating_list;


GaOld::GaOld(Population &population, unsigned dim){
    this->dim = dim;
    this->population = population;
    this->pop_size = population.pop_size;
    mating_list = new double*[2 * pop_size];
    for (unsigned i = 0; i < 2 * pop_size; i++) {
        mating_list[i] = new double[dim];
    }
}

void GaOld::selection_roulette() const {
    // create list of offsets for roulette selection
    auto* offset = new double[pop_size];
    double total_fitness = 0;

    // shift fitness to all positive values and invert it so the minimal value has the highest fitness
    // #pragma omp parallel for default(none) shared(population, pop_size, offset) reduction(+:total_fitness)
    for (unsigned i=0; i<pop_size; i++) {
        // shift fitness to all positive values and invert it so the minimal value has the highest fitness
        population.fitness.fitness_values[i] = population.fitness.max_fitness - population.fitness.fitness_values[i];
        // calculate total fitness
        total_fitness += population.fitness.fitness_values[i];
    }

    offset[0] = population.fitness.fitness_values[0]/total_fitness;
    // calculate offset for roulette selection
    for (unsigned i=1; i<pop_size; i++) {
        offset[i] = offset[i-1] + (population.fitness.fitness_values[i]/total_fitness);
    }

    // do roulette selection
    for (unsigned i=0; i < 2*pop_size; i++) {
        double roulette_random = uniformRealDistribution(gen);
        for (unsigned j=0; j<pop_size; j++) {
            if (roulette_random < offset[j]) {
                for (unsigned k=0; k < dim; k++) {
                    mating_list[i][k] = population.pop[j][k];          // ToDo prevent mating with itself
                }
                break;
            }
        }
    }

    // free memory
    delete[] offset;
}

void GaOld::crossover_uniform() const {
    // #pragma omp parallel for collapse(2) default(none) shared(population, mating_list, pop_size, dim, uniformRealDistribution, gen)
    for (unsigned  i=0; i<pop_size; i++) {
        for (unsigned  j=0; j<dim; j++) {
            if (uniformRealDistribution(gen) < 0.5) {       // choose gene of parent A
                population.pop[i][j] = mating_list[2*i][j];
            } else {                    // choose gene of parent B
                population.pop[i][j] = mating_list[2*i+1][j];
            }
        }
    }
}

void GaOld::evolve(int generations) {
    for (int i=0; i<generations; i++) {
        selection_roulette();
        crossover_uniform();
        // ToDo mutation
        population.fitness.compute_fitness(population.pop);
        // printf("%e \n", population.fitness.min_fitness);
    }
}

void GaOld::clean() const{
    for (unsigned i = 0; i < 2 * pop_size; i++) {
        delete[] mating_list[i];
    }
    delete[] mating_list;
}


