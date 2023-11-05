//
// Created by Marcel Beyer on 22.08.2023.
//
#include <cstdio>
#include "ga.h"
#include <chrono>
#include <limits>
#include <vector>

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

    unsigned percentage = 30;
    runs = 5;
    pop_size = 128;
    unsigned no_sub_pops = 32;  // must be a power of 2

    for (int i = 1; i<=runs ; i++) {
        printf("Run: %i/%i\n", i, runs);
        std::vector<Ga*> ga_vector;

        ga_vector.reserve(no_sub_pops);
        for (int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
            ga_vector.push_back(new Ga(dim, pop_size, -100, 100));
        }

        for(int j = 0; j<4; j++) {
            // first evolution until 50 generations
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
                ga_vector[sub_pop]->evolve(50, false);
            }

            // exchange the best individuals
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop+=2) {
                exchange_percentage(ga_vector[sub_pop]->pop, ga_vector[sub_pop]->fitness, ga_vector[sub_pop+1]->pop, ga_vector[sub_pop+1]->fitness, percentage);
            }

            // generations up to 100
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
                ga_vector[sub_pop]->evolve(50, false);
            }

            // exchange the best individuals
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop+=4) {
                exchange_percentage(ga_vector[sub_pop]->pop, ga_vector[sub_pop]->fitness, ga_vector[sub_pop+2]->pop, ga_vector[sub_pop+2]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+1]->pop, ga_vector[sub_pop+1]->fitness, ga_vector[sub_pop+3]->pop, ga_vector[sub_pop+3]->fitness, percentage);
            }

            // generations up to 150
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
                ga_vector[sub_pop]->evolve(50, false);
            }

            // exchange the best individuals
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop+=8) {
                exchange_percentage(ga_vector[sub_pop]->pop, ga_vector[sub_pop]->fitness, ga_vector[sub_pop+4]->pop, ga_vector[sub_pop+4]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+1]->pop, ga_vector[sub_pop+1]->fitness, ga_vector[sub_pop+5]->pop, ga_vector[sub_pop+5]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+2]->pop, ga_vector[sub_pop+2]->fitness, ga_vector[sub_pop+6]->pop, ga_vector[sub_pop+6]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+3]->pop, ga_vector[sub_pop+3]->fitness, ga_vector[sub_pop+7]->pop, ga_vector[sub_pop+7]->fitness, percentage);
            }

            // generations up to 200
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
                ga_vector[sub_pop]->evolve(50, false);
            }

            // exchange the best individuals
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop+=16) {
                exchange_percentage(ga_vector[sub_pop]->pop, ga_vector[sub_pop]->fitness, ga_vector[sub_pop+8]->pop, ga_vector[sub_pop+8]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+1]->pop, ga_vector[sub_pop+1]->fitness, ga_vector[sub_pop+9]->pop, ga_vector[sub_pop+9]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+2]->pop, ga_vector[sub_pop+2]->fitness, ga_vector[sub_pop+10]->pop, ga_vector[sub_pop+10]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+3]->pop, ga_vector[sub_pop+3]->fitness, ga_vector[sub_pop+11]->pop, ga_vector[sub_pop+11]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+4]->pop, ga_vector[sub_pop+4]->fitness, ga_vector[sub_pop+12]->pop, ga_vector[sub_pop+12]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+5]->pop, ga_vector[sub_pop+5]->fitness, ga_vector[sub_pop+13]->pop, ga_vector[sub_pop+13]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+6]->pop, ga_vector[sub_pop+6]->fitness, ga_vector[sub_pop+14]->pop, ga_vector[sub_pop+14]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+7]->pop, ga_vector[sub_pop+7]->fitness, ga_vector[sub_pop+15]->pop, ga_vector[sub_pop+15]->fitness, percentage);
            }

            // generations up to 250
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
                ga_vector[sub_pop]->evolve(50, false);
            }

            // exchange the best individuals
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop+=32) {
                exchange_percentage(ga_vector[sub_pop]->pop, ga_vector[sub_pop]->fitness, ga_vector[sub_pop+16]->pop, ga_vector[sub_pop+16]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+1]->pop, ga_vector[sub_pop+1]->fitness, ga_vector[sub_pop+17]->pop, ga_vector[sub_pop+17]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+2]->pop, ga_vector[sub_pop+2]->fitness, ga_vector[sub_pop+18]->pop, ga_vector[sub_pop+18]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+3]->pop, ga_vector[sub_pop+3]->fitness, ga_vector[sub_pop+19]->pop, ga_vector[sub_pop+19]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+4]->pop, ga_vector[sub_pop+4]->fitness, ga_vector[sub_pop+20]->pop, ga_vector[sub_pop+20]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+5]->pop, ga_vector[sub_pop+5]->fitness, ga_vector[sub_pop+21]->pop, ga_vector[sub_pop+21]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+6]->pop, ga_vector[sub_pop+6]->fitness, ga_vector[sub_pop+22]->pop, ga_vector[sub_pop+22]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+7]->pop, ga_vector[sub_pop+7]->fitness, ga_vector[sub_pop+23]->pop, ga_vector[sub_pop+23]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+8]->pop, ga_vector[sub_pop+8]->fitness, ga_vector[sub_pop+24]->pop, ga_vector[sub_pop+24]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+9]->pop, ga_vector[sub_pop+9]->fitness, ga_vector[sub_pop+25]->pop, ga_vector[sub_pop+25]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+10]->pop, ga_vector[sub_pop+10]->fitness, ga_vector[sub_pop+26]->pop, ga_vector[sub_pop+26]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+11]->pop, ga_vector[sub_pop+11]->fitness, ga_vector[sub_pop+27]->pop, ga_vector[sub_pop+27]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+12]->pop, ga_vector[sub_pop+12]->fitness, ga_vector[sub_pop+28]->pop, ga_vector[sub_pop+28]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+13]->pop, ga_vector[sub_pop+13]->fitness, ga_vector[sub_pop+29]->pop, ga_vector[sub_pop+29]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+14]->pop, ga_vector[sub_pop+14]->fitness, ga_vector[sub_pop+30]->pop, ga_vector[sub_pop+30]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+15]->pop, ga_vector[sub_pop+15]->fitness, ga_vector[sub_pop+31]->pop, ga_vector[sub_pop+31]->fitness, percentage);
            }
        }

        for (int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
            ga_vector[sub_pop]->~Ga();
        }
    }


# endif

#if 0
    runs = 5;
    pop_size = 256;
    for(unsigned i = 1; i <= runs; i++) {
        printf("Run: %i/%i at %i individuals\n", i, runs, pop_size);
        Ga* ga0 = new Ga(dim, pop_size, -100, 100);
        ga0->evolve(1000, false);

        ga0->~Ga();
    }

    runs = 5;
    pop_size = 128;
    for(unsigned i = 1; i <= runs; i++) {
        printf("Run: %i/%i at %i individuals\n", i, runs, pop_size);
        Ga* ga0 = new Ga(dim, pop_size, -100, 100);
        ga0->evolve(1000, false);

        ga0->~Ga();
    }
#endif

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    printf("Time taken by function: %lld seconds", duration.count());
    return 0;
}