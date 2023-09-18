//
// Created by Marcel Beyer on 22.08.2023.
//
#include <cstdio>
#include "ga.h"
#include <chrono>
#include <cstdlib>
#include <limits>

using namespace std::chrono;

unsigned dim = 1000;
unsigned pop_size = 5000;

unsigned min_index(const double* array, unsigned size) {
    unsigned min_index = 0;
    for (unsigned i = 1; i < size; i++) {
        if (array[i] < array[min_index]) {
            min_index = i;
        }
    }
    return min_index;
}

void copy_pop_element(double** source, unsigned source_index, double** target, unsigned target_index) {
    for (unsigned i = 0; i < dim; i++) {
        target[target_index][i] = source[source_index][i];
    }
}

void copy_best_individuals(double** a_pop, double* a_fit, double** b_pop, double* b_fit, double** target) {
    for (unsigned i = 0; i < pop_size; i++) {
        unsigned min_a = min_index(a_fit, pop_size);
        unsigned min_b = min_index(b_fit, pop_size);
        if(a_fit[min_a] < b_fit[min_b]) {
            copy_pop_element(a_pop, min_a, target, i);
            a_fit[min_a] = std::numeric_limits<double>::max();
        } else {
            copy_pop_element(b_pop, min_b, target, i);
            b_fit[min_b] = std::numeric_limits<double>::max();
        }
    }
}

int main() {
    auto start = high_resolution_clock::now();

    Ga* ga0 = new Ga(dim, pop_size, -100, 100);
    ga0->evolve(1000, true);
    // Ga* ga1 = new Ga(dim, pop_size, -100, 100);

    /*
    for(unsigned k = 0; k < 3; k++) {
        ga0->evolve(350, true);
        ga1->evolve(350, true);

        auto* fitness_copy = new double[pop_size];
        for(unsigned i = 0; i < pop_size; i++) {
            fitness_copy[i] = ga0->fitness[i];
        }

        auto** best_individuals = new double*[pop_size];
        for(unsigned i = 0; i < pop_size; i++) {
            best_individuals[i] = new double[dim];
        }

        copy_best_individuals(ga0->pop, ga0->fitness, ga1->pop, ga1->fitness, best_individuals);

        for(unsigned i = 0; i < pop_size; i++) {
            copy_pop_element(best_individuals, i, ga0->pop, i);
            copy_pop_element(best_individuals, i, ga1->pop, i);
        }

        // free allocated memory
        for(unsigned i = 0; i < pop_size; i++) {
            delete[] best_individuals[i];
        }
        delete[] best_individuals;
        delete[] fitness_copy;
    }
    */

    // free allocated memory
    ga0->clean();
    // ga1->clean();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    printf("Time taken by function: %lld seconds", duration.count());
    return 0;
}