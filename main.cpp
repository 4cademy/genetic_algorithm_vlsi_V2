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
    for (int i = 1; i<=runs ; i++) {
        printf("Run: %i/%i\n", i, runs);
        Ga *ga0 = new Ga(dim, pop_size, -100, 100);
        Ga *ga1 = new Ga(dim, pop_size, -100, 100);
        Ga *ga2 = new Ga(dim, pop_size, -100, 100);
        Ga *ga3 = new Ga(dim, pop_size, -100, 100);
        Ga *ga4 = new Ga(dim, pop_size, -100, 100);
        Ga *ga5 = new Ga(dim, pop_size, -100, 100);
        Ga *ga6 = new Ga(dim, pop_size, -100, 100);
        Ga *ga7 = new Ga(dim, pop_size, -100, 100);
        Ga *ga8 = new Ga(dim, pop_size, -100, 100);
        Ga *ga9 = new Ga(dim, pop_size, -100, 100);
        Ga *ga10 = new Ga(dim, pop_size, -100, 100);
        Ga *ga11 = new Ga(dim, pop_size, -100, 100);
        Ga *ga12 = new Ga(dim, pop_size, -100, 100);
        Ga *ga13 = new Ga(dim, pop_size, -100, 100);
        Ga *ga14 = new Ga(dim, pop_size, -100, 100);
        Ga *ga15 = new Ga(dim, pop_size, -100, 100);
        Ga *ga16 = new Ga(dim, pop_size, -100, 100);
        Ga *ga17 = new Ga(dim, pop_size, -100, 100);
        Ga *ga18 = new Ga(dim, pop_size, -100, 100);
        Ga *ga19 = new Ga(dim, pop_size, -100, 100);
        Ga *ga20 = new Ga(dim, pop_size, -100, 100);
        Ga *ga21 = new Ga(dim, pop_size, -100, 100);
        Ga *ga22 = new Ga(dim, pop_size, -100, 100);
        Ga *ga23 = new Ga(dim, pop_size, -100, 100);
        Ga *ga24 = new Ga(dim, pop_size, -100, 100);
        Ga *ga25 = new Ga(dim, pop_size, -100, 100);
        Ga *ga26 = new Ga(dim, pop_size, -100, 100);
        Ga *ga27 = new Ga(dim, pop_size, -100, 100);
        Ga *ga28 = new Ga(dim, pop_size, -100, 100);
        Ga *ga29 = new Ga(dim, pop_size, -100, 100);
        Ga *ga30 = new Ga(dim, pop_size, -100, 100);
        Ga *ga31 = new Ga(dim, pop_size, -100, 100);

        for(int j = 0; j<4; j++) {
            // first evolution until 50 generations
            ga0->evolve(50, false);
            ga1->evolve(50, false);
            ga2->evolve(50, false);
            ga3->evolve(50, false);
            ga4->evolve(50, false);
            ga5->evolve(50, false);
            ga6->evolve(50, false);
            ga7->evolve(50, false);
            ga8->evolve(50, false);
            ga9->evolve(50, false);
            ga10->evolve(50, false);
            ga11->evolve(50, false);
            ga12->evolve(50, false);
            ga13->evolve(50, false);
            ga14->evolve(50, false);
            ga15->evolve(50, false);
            ga16->evolve(50, false);
            ga17->evolve(50, false);
            ga18->evolve(50, false);
            ga19->evolve(50, false);
            ga20->evolve(50, false);
            ga21->evolve(50, false);
            ga22->evolve(50, false);
            ga23->evolve(50, false);
            ga24->evolve(50, false);
            ga25->evolve(50, false);
            ga26->evolve(50, false);
            ga27->evolve(50, false);
            ga28->evolve(50, false);
            ga29->evolve(50, false);
            ga30->evolve(50, false);
            ga31->evolve(50, false);

            exchange_percentage(ga0->pop, ga0->fitness, ga1->pop, ga1->fitness, percentage);

            exchange_percentage(ga2->pop, ga2->fitness, ga3->pop, ga3->fitness, percentage);

            exchange_percentage(ga4->pop, ga4->fitness, ga5->pop, ga5->fitness, percentage);

            exchange_percentage(ga6->pop, ga6->fitness, ga7->pop, ga7->fitness, percentage);

            exchange_percentage(ga8->pop, ga8->fitness, ga9->pop, ga9->fitness, percentage);

            exchange_percentage(ga10->pop, ga10->fitness, ga11->pop, ga11->fitness, percentage);

            exchange_percentage(ga12->pop, ga12->fitness, ga13->pop, ga13->fitness, percentage);

            exchange_percentage(ga14->pop, ga14->fitness, ga15->pop, ga15->fitness, percentage);

            exchange_percentage(ga16->pop, ga16->fitness, ga17->pop, ga17->fitness, percentage);

            exchange_percentage(ga18->pop, ga18->fitness, ga19->pop, ga19->fitness, percentage);

            exchange_percentage(ga20->pop, ga20->fitness, ga21->pop, ga21->fitness, percentage);

            exchange_percentage(ga22->pop, ga22->fitness, ga23->pop, ga23->fitness, percentage);

            exchange_percentage(ga24->pop, ga24->fitness, ga25->pop, ga25->fitness, percentage);

            exchange_percentage(ga26->pop, ga26->fitness, ga27->pop, ga27->fitness, percentage);

            exchange_percentage(ga28->pop, ga28->fitness, ga29->pop, ga29->fitness, percentage);

            exchange_percentage(ga30->pop, ga30->fitness, ga31->pop, ga31->fitness, percentage);

            // generations up to 100
            ga0->evolve(50, false);
            ga1->evolve(50, false);
            ga2->evolve(50, false);
            ga3->evolve(50, false);
            ga4->evolve(50, false);
            ga5->evolve(50, false);
            ga6->evolve(50, false);
            ga7->evolve(50, false);
            ga8->evolve(50, false);
            ga9->evolve(50, false);
            ga10->evolve(50, false);
            ga11->evolve(50, false);
            ga12->evolve(50, false);
            ga13->evolve(50, false);
            ga14->evolve(50, false);
            ga15->evolve(50, false);
            ga16->evolve(50, false);
            ga17->evolve(50, false);
            ga18->evolve(50, false);
            ga19->evolve(50, false);
            ga20->evolve(50, false);
            ga21->evolve(50, false);
            ga22->evolve(50, false);
            ga23->evolve(50, false);
            ga24->evolve(50, false);
            ga25->evolve(50, false);
            ga26->evolve(50, false);
            ga27->evolve(50, false);
            ga28->evolve(50, false);
            ga29->evolve(50, false);
            ga30->evolve(50, false);
            ga31->evolve(50, false);

            exchange_percentage(ga0->pop, ga0->fitness, ga2->pop, ga2->fitness, percentage);
            exchange_percentage(ga1->pop, ga1->fitness, ga3->pop, ga3->fitness, percentage);

            exchange_percentage(ga4->pop, ga4->fitness, ga6->pop, ga6->fitness, percentage);
            exchange_percentage(ga5->pop, ga5->fitness, ga7->pop, ga7->fitness, percentage);

            exchange_percentage(ga8->pop, ga8->fitness, ga10->pop, ga10->fitness, percentage);
            exchange_percentage(ga9->pop, ga9->fitness, ga11->pop, ga11->fitness, percentage);

            exchange_percentage(ga12->pop, ga12->fitness, ga14->pop, ga14->fitness, percentage);
            exchange_percentage(ga13->pop, ga13->fitness, ga15->pop, ga15->fitness, percentage);

            exchange_percentage(ga16->pop, ga16->fitness, ga18->pop, ga18->fitness, percentage);
            exchange_percentage(ga17->pop, ga17->fitness, ga19->pop, ga19->fitness, percentage);

            exchange_percentage(ga20->pop, ga20->fitness, ga22->pop, ga22->fitness, percentage);
            exchange_percentage(ga21->pop, ga21->fitness, ga23->pop, ga23->fitness, percentage);

            exchange_percentage(ga24->pop, ga24->fitness, ga26->pop, ga26->fitness, percentage);
            exchange_percentage(ga25->pop, ga25->fitness, ga27->pop, ga27->fitness, percentage);

            exchange_percentage(ga28->pop, ga28->fitness, ga30->pop, ga30->fitness, percentage);
            exchange_percentage(ga29->pop, ga29->fitness, ga31->pop, ga31->fitness, percentage);

            // generations up to 150
            ga0->evolve(50, false);
            ga1->evolve(50, false);
            ga2->evolve(50, false);
            ga3->evolve(50, false);
            ga4->evolve(50, false);
            ga5->evolve(50, false);
            ga6->evolve(50, false);
            ga7->evolve(50, false);
            ga8->evolve(50, false);
            ga9->evolve(50, false);
            ga10->evolve(50, false);
            ga11->evolve(50, false);
            ga12->evolve(50, false);
            ga13->evolve(50, false);
            ga14->evolve(50, false);
            ga15->evolve(50, false);
            ga16->evolve(50, false);
            ga17->evolve(50, false);
            ga18->evolve(50, false);
            ga19->evolve(50, false);
            ga20->evolve(50, false);
            ga21->evolve(50, false);
            ga22->evolve(50, false);
            ga23->evolve(50, false);
            ga24->evolve(50, false);
            ga25->evolve(50, false);
            ga26->evolve(50, false);
            ga27->evolve(50, false);
            ga28->evolve(50, false);
            ga29->evolve(50, false);
            ga30->evolve(50, false);
            ga31->evolve(50, false);

            exchange_percentage(ga0->pop, ga0->fitness, ga4->pop, ga4->fitness, percentage);
            exchange_percentage(ga1->pop, ga1->fitness, ga5->pop, ga5->fitness, percentage);
            exchange_percentage(ga2->pop, ga2->fitness, ga6->pop, ga6->fitness, percentage);
            exchange_percentage(ga3->pop, ga3->fitness, ga7->pop, ga7->fitness, percentage);

            exchange_percentage(ga8->pop, ga8->fitness, ga12->pop, ga12->fitness, percentage);
            exchange_percentage(ga9->pop, ga9->fitness, ga13->pop, ga13->fitness, percentage);
            exchange_percentage(ga10->pop, ga10->fitness, ga14->pop, ga14->fitness, percentage);
            exchange_percentage(ga11->pop, ga11->fitness, ga15->pop, ga15->fitness, percentage);

            exchange_percentage(ga16->pop, ga16->fitness, ga20->pop, ga20->fitness, percentage);
            exchange_percentage(ga17->pop, ga17->fitness, ga21->pop, ga21->fitness, percentage);
            exchange_percentage(ga18->pop, ga18->fitness, ga22->pop, ga22->fitness, percentage);
            exchange_percentage(ga19->pop, ga19->fitness, ga23->pop, ga23->fitness, percentage);

            exchange_percentage(ga24->pop, ga24->fitness, ga28->pop, ga28->fitness, percentage);
            exchange_percentage(ga25->pop, ga25->fitness, ga29->pop, ga29->fitness, percentage);
            exchange_percentage(ga26->pop, ga26->fitness, ga30->pop, ga30->fitness, percentage);
            exchange_percentage(ga27->pop, ga27->fitness, ga31->pop, ga31->fitness, percentage);

            // generations up to 200
            ga0->evolve(50, false);
            ga1->evolve(50, false);
            ga2->evolve(50, false);
            ga3->evolve(50, false);
            ga4->evolve(50, false);
            ga5->evolve(50, false);
            ga6->evolve(50, false);
            ga7->evolve(50, false);
            ga8->evolve(50, false);
            ga9->evolve(50, false);
            ga10->evolve(50, false);
            ga11->evolve(50, false);
            ga12->evolve(50, false);
            ga13->evolve(50, false);
            ga14->evolve(50, false);
            ga15->evolve(50, false);
            ga16->evolve(50, false);
            ga17->evolve(50, false);
            ga18->evolve(50, false);
            ga19->evolve(50, false);
            ga20->evolve(50, false);
            ga21->evolve(50, false);
            ga22->evolve(50, false);
            ga23->evolve(50, false);
            ga24->evolve(50, false);
            ga25->evolve(50, false);
            ga26->evolve(50, false);
            ga27->evolve(50, false);
            ga28->evolve(50, false);
            ga29->evolve(50, false);
            ga30->evolve(50, false);
            ga31->evolve(50, false);

            exchange_percentage(ga0->pop, ga0->fitness, ga8->pop, ga8->fitness, percentage);
            exchange_percentage(ga1->pop, ga1->fitness, ga9->pop, ga9->fitness, percentage);
            exchange_percentage(ga2->pop, ga2->fitness, ga10->pop, ga10->fitness, percentage);
            exchange_percentage(ga3->pop, ga3->fitness, ga11->pop, ga11->fitness, percentage);
            exchange_percentage(ga4->pop, ga4->fitness, ga12->pop, ga12->fitness, percentage);
            exchange_percentage(ga5->pop, ga5->fitness, ga13->pop, ga13->fitness, percentage);
            exchange_percentage(ga6->pop, ga6->fitness, ga14->pop, ga14->fitness, percentage);
            exchange_percentage(ga7->pop, ga7->fitness, ga15->pop, ga15->fitness, percentage);

            exchange_percentage(ga16->pop, ga16->fitness, ga24->pop, ga24->fitness, percentage);
            exchange_percentage(ga17->pop, ga17->fitness, ga25->pop, ga25->fitness, percentage);
            exchange_percentage(ga18->pop, ga18->fitness, ga26->pop, ga26->fitness, percentage);
            exchange_percentage(ga19->pop, ga19->fitness, ga27->pop, ga27->fitness, percentage);
            exchange_percentage(ga20->pop, ga20->fitness, ga28->pop, ga28->fitness, percentage);
            exchange_percentage(ga21->pop, ga21->fitness, ga29->pop, ga29->fitness, percentage);
            exchange_percentage(ga22->pop, ga22->fitness, ga30->pop, ga30->fitness, percentage);
            exchange_percentage(ga23->pop, ga23->fitness, ga31->pop, ga31->fitness, percentage);

            // generations up to 250
            ga0->evolve(50, false);
            ga1->evolve(50, false);
            ga2->evolve(50, false);
            ga3->evolve(50, false);
            ga4->evolve(50, false);
            ga5->evolve(50, false);
            ga6->evolve(50, false);
            ga7->evolve(50, false);
            ga8->evolve(50, false);
            ga9->evolve(50, false);
            ga10->evolve(50, false);
            ga11->evolve(50, false);
            ga12->evolve(50, false);
            ga13->evolve(50, false);
            ga14->evolve(50, false);
            ga15->evolve(50, false);
            ga16->evolve(50, false);
            ga17->evolve(50, false);
            ga18->evolve(50, false);
            ga19->evolve(50, false);
            ga20->evolve(50, false);
            ga21->evolve(50, false);
            ga22->evolve(50, false);
            ga23->evolve(50, false);
            ga24->evolve(50, false);
            ga25->evolve(50, false);
            ga26->evolve(50, false);
            ga27->evolve(50, false);
            ga28->evolve(50, false);
            ga29->evolve(50, false);
            ga30->evolve(50, false);
            ga31->evolve(50, false);

            exchange_percentage(ga0->pop, ga0->fitness, ga16->pop, ga16->fitness, percentage);
            exchange_percentage(ga1->pop, ga1->fitness, ga17->pop, ga17->fitness, percentage);
            exchange_percentage(ga2->pop, ga2->fitness, ga18->pop, ga18->fitness, percentage);
            exchange_percentage(ga3->pop, ga3->fitness, ga19->pop, ga19->fitness, percentage);
            exchange_percentage(ga4->pop, ga4->fitness, ga20->pop, ga20->fitness, percentage);
            exchange_percentage(ga5->pop, ga5->fitness, ga21->pop, ga21->fitness, percentage);
            exchange_percentage(ga6->pop, ga6->fitness, ga22->pop, ga22->fitness, percentage);
            exchange_percentage(ga7->pop, ga7->fitness, ga23->pop, ga23->fitness, percentage);
            exchange_percentage(ga8->pop, ga8->fitness, ga24->pop, ga24->fitness, percentage);
            exchange_percentage(ga9->pop, ga9->fitness, ga25->pop, ga25->fitness, percentage);
            exchange_percentage(ga10->pop, ga10->fitness, ga26->pop, ga26->fitness, percentage);
            exchange_percentage(ga11->pop, ga11->fitness, ga27->pop, ga27->fitness, percentage);
            exchange_percentage(ga12->pop, ga12->fitness, ga28->pop, ga28->fitness, percentage);
            exchange_percentage(ga13->pop, ga13->fitness, ga29->pop, ga29->fitness, percentage);
            exchange_percentage(ga14->pop, ga14->fitness, ga30->pop, ga30->fitness, percentage);
            exchange_percentage(ga15->pop, ga15->fitness, ga31->pop, ga31->fitness, percentage);
        }

        ga0->~Ga();
        ga1->~Ga();
        ga2->~Ga();
        ga3->~Ga();
        ga4->~Ga();
        ga5->~Ga();
        ga6->~Ga();
        ga7->~Ga();
        ga8->~Ga();
        ga9->~Ga();
        ga10->~Ga();
        ga11->~Ga();
        ga12->~Ga();
        ga13->~Ga();
        ga14->~Ga();
        ga15->~Ga();
        ga16->~Ga();
        ga17->~Ga();
        ga18->~Ga();
        ga19->~Ga();
        ga20->~Ga();
        ga21->~Ga();
        ga22->~Ga();
        ga23->~Ga();
        ga24->~Ga();
        ga25->~Ga();
        ga26->~Ga();
        ga27->~Ga();
        ga28->~Ga();
        ga29->~Ga();
        ga30->~Ga();
        ga31->~Ga();
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