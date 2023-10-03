//
// Created by Marcel Beyer on 22.08.2023.
//
#include <cstdio>
#include "ga.h"
#include <chrono>
#include <limits>

using namespace std::chrono;

unsigned dim = 1'000;
unsigned pop_size = 10'000;
unsigned runs = 0;

enum minmax {
    MIN,
    MAX
};

/*
 * this function returns the index of the min/max fitness value
 * 'min_or_max' specifies whether the min or max index should be returned
 */
unsigned min_max_index(const float* fitness, unsigned size, minmax min_or_max) {
    unsigned index = 0;
    for (unsigned i = 1; i < size; i++) {
        if (min_or_max == MAX) {
            if (fitness[i] > fitness[index]) {
                index = i;
            }
        } else{
            if (fitness[i] < fitness[index]) {
                index = i;
            }
        }
    }
    return index;
}

/*
 * this function returns a dynamic list of the indices of the min/max fitness values
 * 'amount' specifies the number of indices to be returned
 * 'min_or_max' specifies whether the min or max indices should be returned
 */
unsigned* min_max_index_list(const float* fitness, unsigned size, unsigned amount, minmax min_or_max) {
    auto* index_list = new unsigned [amount];

    // make temp copy of fitness
    auto* fitness_copy = new float[size];
    for (unsigned i = 0; i < size; i++) {
        fitness_copy[i] = fitness[i];
    }

    if(min_or_max == MAX) {
        for (unsigned i = 0; i < amount; i++) {
            index_list[i] = min_max_index(fitness_copy, size, MAX);
            fitness_copy[index_list[i]] = std::numeric_limits<float>::min();
        }
    } else {
        for (unsigned i = 0; i < amount; i++) {
            index_list[i] = min_max_index(fitness_copy, size, MIN);
            fitness_copy[index_list[i]] = std::numeric_limits<float>::max();
        }
    }

    delete[] fitness_copy;
    return index_list;
}

/*
 * this function copies one individual from 'source' at 'source_index' to 'target' at 'target_index'
 */
void copy_individual(float** source, unsigned source_index, float** target, unsigned target_index) {
    for (unsigned i = 0; i < dim; i++) {
        target[target_index][i] = source[source_index][i];
    }
}

// this function creates a new population by copying the best individuals from a_pop and b_pop
void copy_best_individuals(float** a_pop, const float* a_fit, float** b_pop, const float* b_fit, float** target) {
    // make temp copies of a_fit and b_fit
    auto* a_fit_copy = new float[pop_size];
    for (unsigned i = 0; i < pop_size; i++) {
        a_fit_copy[i] = a_fit[i];
    }
    auto* b_fit_copy = new float[pop_size];
    for (unsigned i = 0; i < pop_size; i++) {
        b_fit_copy[i] = b_fit[i];
    }

    for (unsigned i = 0; i < pop_size; i++) {
        unsigned min_a = min_max_index(a_fit_copy, pop_size, MIN);
        unsigned min_b = min_max_index(b_fit_copy, pop_size, MIN);
        if(a_fit_copy[min_a] < b_fit_copy[min_b]) {
            copy_individual(a_pop, min_a, target, i);
            a_fit_copy[min_a] = std::numeric_limits<float>::max();
        } else {
            copy_individual(b_pop, min_b, target, i);
            b_fit_copy[min_b] = std::numeric_limits<float>::max();
        }
    }

    delete[] a_fit_copy;
    delete[] b_fit_copy;
}

// this function creates a new population by exchanging a percentage of the best individuals from a_pop and b_pop
void exchange_percentage(float** a_pop, float* a_fit, float** b_pop, float* b_fit, unsigned percentage) {
    unsigned exchange_amount = pop_size * percentage / 100;
    unsigned* a_min_list = min_max_index_list(a_fit, pop_size, exchange_amount, MIN);
    unsigned* a_max_list = min_max_index_list(a_fit, pop_size, exchange_amount, MAX);
    unsigned* b_min_list = min_max_index_list(b_fit, pop_size, exchange_amount, MIN);
    unsigned* b_max_list = min_max_index_list(b_fit, pop_size, exchange_amount, MAX);

    for (unsigned i = 0; i < exchange_amount; i++) {
        copy_individual(a_pop, a_min_list[i], b_pop, b_max_list[i]);
        copy_individual(b_pop, b_min_list[i], a_pop, a_max_list[i]);
    }
}

int main() {
    auto start = high_resolution_clock::now();

#if 0
    runs = 0;
    pop_size = 10'000;
    for(unsigned i = 1; i <= runs; i++) {
        printf("Run: %i/%i at %i individuals\n", i, runs, pop_size);
        Ga* ga0 = new Ga(dim, pop_size, -100, 100);
        ga0->evolve(1000, false);
        ga0->~Ga();
    }
    runs = 0;
    pop_size = 5'000;
    for(unsigned i = 1; i <= runs; i++) {
        printf("Run: %i/%i at %i individuals\n", i, runs, pop_size);
        Ga* ga0 = new Ga(dim, pop_size, -100, 100);
        ga0->evolve(1000, false);
        ga0->~Ga();
    }
    runs = 3;
    pop_size = 20'000;
    for(unsigned i = 1; i <= runs; i++) {
        printf("Run: %i/%i at %i individuals\n", i, runs, pop_size);
        Ga* ga0 = new Ga(dim, pop_size, -100, 100);
        ga0->evolve(1000, false);
        ga0->~Ga();
    }
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
            copy_individual(best_individuals, i, ga0->pop, i);
            copy_individual(best_individuals, i, ga1->pop, i);
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

# if 1
    unsigned percentage = 10;
    runs = 2;
    pop_size = 5'000;
    printf("Preceding %i percent exchange:\n", percentage);
    for (int i = 1; i<=runs ; i++) {
        printf("Run: %i/%i\n", i, runs);
        Ga *ga0 = new Ga(dim, pop_size, -100, 100);
        Ga *ga1 = new Ga(dim, pop_size, -100, 100);

        for (unsigned k = 0; k < 1; k++) {
            ga0->evolve(600, false);
            ga1->evolve(600, false);

            exchange_percentage(ga0->pop, ga0->fitness, ga1->pop, ga1->fitness, percentage);
        }
        ga0->evolve(400, false);
        ga1->evolve(400, false);

        ga0->~Ga();
        ga1->~Ga();
    }

    percentage = 20;
    runs = 2;
    printf("Preceding %i percent exchange:\n", percentage);
    for (int i = 1; i<=runs ; i++) {
        printf("Run: %i/%i\n", i, runs);
        Ga *ga0 = new Ga(dim, pop_size, -100, 100);
        Ga *ga1 = new Ga(dim, pop_size, -100, 100);

        for (unsigned k = 0; k < 1; k++) {
            ga0->evolve(600, false);
            ga1->evolve(600, false);

            exchange_percentage(ga0->pop, ga0->fitness, ga1->pop, ga1->fitness, percentage);
        }
        ga0->evolve(400, false);
        ga1->evolve(400, false);

        ga0->~Ga();
        ga1->~Ga();
    }

    percentage = 30;
    runs = 2;
    printf("Preceding %i percent exchange:\n", percentage);
    for (int i = 1; i<=runs ; i++) {
        printf("Run: %i/%i\n", i, runs);
        Ga *ga0 = new Ga(dim, pop_size, -100, 100);
        Ga *ga1 = new Ga(dim, pop_size, -100, 100);

        for (unsigned k = 0; k < 1; k++) {
            ga0->evolve(600, false);
            ga1->evolve(600, false);

            exchange_percentage(ga0->pop, ga0->fitness, ga1->pop, ga1->fitness, percentage);
        }
        ga0->evolve(400, false);
        ga1->evolve(400, false);

        ga0->~Ga();
        ga1->~Ga();
    }
# endif

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    printf("Time taken by function: %lld seconds", duration.count());
    return 0;
}