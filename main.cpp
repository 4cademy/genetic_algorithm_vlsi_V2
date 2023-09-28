//
// Created by Marcel Beyer on 22.08.2023.
//
#include <cstdio>
#include "ga.h"
#include <chrono>
#include <cstdlib>
#include <limits>

using namespace std::chrono;

unsigned dim = 1'000;
unsigned pop_size = 5'000;

unsigned min_index(const float* array, unsigned size) {
    unsigned min_index = 0;
    for (unsigned i = 1; i < size; i++) {
        if (array[i] < array[min_index]) {
            min_index = i;
        }
    }
    return min_index;
}

void copy_pop_element(float** source, unsigned source_index, float** target, unsigned target_index) {
    for (unsigned i = 0; i < dim; i++) {
        target[target_index][i] = source[source_index][i];
    }
}

void copy_best_individuals(float** a_pop, float* a_fit, float** b_pop, float* b_fit, float** target) {
    for (unsigned i = 0; i < pop_size; i++) {
        unsigned min_a = min_index(a_fit, pop_size);
        unsigned min_b = min_index(b_fit, pop_size);
        if(a_fit[min_a] < b_fit[min_b]) {
            copy_pop_element(a_pop, min_a, target, i);
            a_fit[min_a] = std::numeric_limits<float>::max();
        } else {
            copy_pop_element(b_pop, min_b, target, i);
            b_fit[min_b] = std::numeric_limits<float>::max();
        }
    }
}

int main() {
    auto start = high_resolution_clock::now();

#if 1
    Ga* ga0 = new Ga(dim, pop_size, -100, 100);
    ga0->evolve(1000, false);
    ga0->~Ga();
#endif

# if 0
    Ga* ga0 = new Ga(dim, pop_size, -100, 100);
    Ga* ga1 = new Ga(dim, pop_size, -100, 100);

    for(unsigned k = 0; k < 1; k++) {
        ga0->evolve(600, false);
        ga1->evolve(600, false);

        auto* fitness_copy = new float[pop_size];
        for(unsigned i = 0; i < pop_size; i++) {
            fitness_copy[i] = ga0->fitness[i];
        }

        auto** best_individuals = new float*[pop_size];
        for(unsigned i = 0; i < pop_size; i++) {
            best_individuals[i] = new float[dim];
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
    ga0->evolve(400, false);
    ga1->evolve(400, false);

    ga0->~Ga();
    ga1->~Ga();
# endif

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    printf("Time taken by function: %lld seconds", duration.count());
    return 0;
}