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
unsigned runs = 5;
float min_gene = -100;
float max_gene = 100;

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

# if 1
    runs = 5;
    pop_size = 4096;
    for (int i = 1; i<=runs ; i++) {
        printf("Run: %i/%i\n", i, runs);
        Ga *ga0 = new Ga(dim, pop_size, min_gene, max_gene);

        ga0->mutation_rate = 0.0001/100;  // 0.0001% mutation rate
        ga0->mutation_deviation = (max_gene-min_gene)/12;  // (max_gene-min_gene)/6 -> 99.8% of the values are in the range [min_gene, max_gene]
        ga0->evolve(1000, false);

        ga0->~Ga();
    }

# endif

# if 0

    unsigned percentage = 30;
    runs = 5;
    pop_size = 128;
    unsigned no_sub_pops = 128;  // must be a power of 2

    for (int i = 1; i<=runs ; i++) {
        printf("Run: %i/%i\n", i, runs);
        std::vector<Ga*> ga_vector;

        ga_vector.reserve(no_sub_pops);
        for (int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
            ga_vector.push_back(new Ga(dim, pop_size, -100, 100));
        }

        for(int j = 0; j<5; j++) {
            // first evolution until 30 generations
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
                ga_vector[sub_pop]->evolve(30, false);
            }

            // exchange the best individuals
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop+=2) {
                exchange_percentage(ga_vector[sub_pop]->pop, ga_vector[sub_pop]->fitness, ga_vector[sub_pop+1]->pop, ga_vector[sub_pop+1]->fitness, percentage);
            }

            // generations up to 60
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
                ga_vector[sub_pop]->evolve(30, false);
            }

            // exchange the best individuals
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop+=4) {
                exchange_percentage(ga_vector[sub_pop]->pop, ga_vector[sub_pop]->fitness, ga_vector[sub_pop+2]->pop, ga_vector[sub_pop+2]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+1]->pop, ga_vector[sub_pop+1]->fitness, ga_vector[sub_pop+3]->pop, ga_vector[sub_pop+3]->fitness, percentage);
            }

            // generations up to 90
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
                ga_vector[sub_pop]->evolve(30, false);
            }

            // exchange the best individuals
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop+=8) {
                exchange_percentage(ga_vector[sub_pop]->pop, ga_vector[sub_pop]->fitness, ga_vector[sub_pop+4]->pop, ga_vector[sub_pop+4]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+1]->pop, ga_vector[sub_pop+1]->fitness, ga_vector[sub_pop+5]->pop, ga_vector[sub_pop+5]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+2]->pop, ga_vector[sub_pop+2]->fitness, ga_vector[sub_pop+6]->pop, ga_vector[sub_pop+6]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+3]->pop, ga_vector[sub_pop+3]->fitness, ga_vector[sub_pop+7]->pop, ga_vector[sub_pop+7]->fitness, percentage);
            }

            // generations up to 120
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
                ga_vector[sub_pop]->evolve(30, false);
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

            // generations up to 150
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
                ga_vector[sub_pop]->evolve(30, false);
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

            // generations up to 180
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
                ga_vector[sub_pop]->evolve(30, false);
            }

            // exchange the best individuals
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop+=64) {
                exchange_percentage(ga_vector[sub_pop]->pop, ga_vector[sub_pop]->fitness, ga_vector[sub_pop+32]->pop, ga_vector[sub_pop+32]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+1]->pop, ga_vector[sub_pop+1]->fitness, ga_vector[sub_pop+33]->pop, ga_vector[sub_pop+33]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+2]->pop, ga_vector[sub_pop+2]->fitness, ga_vector[sub_pop+34]->pop, ga_vector[sub_pop+34]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+3]->pop, ga_vector[sub_pop+3]->fitness, ga_vector[sub_pop+35]->pop, ga_vector[sub_pop+35]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+4]->pop, ga_vector[sub_pop+4]->fitness, ga_vector[sub_pop+36]->pop, ga_vector[sub_pop+36]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+5]->pop, ga_vector[sub_pop+5]->fitness, ga_vector[sub_pop+37]->pop, ga_vector[sub_pop+37]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+6]->pop, ga_vector[sub_pop+6]->fitness, ga_vector[sub_pop+38]->pop, ga_vector[sub_pop+38]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+7]->pop, ga_vector[sub_pop+7]->fitness, ga_vector[sub_pop+39]->pop, ga_vector[sub_pop+39]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+8]->pop, ga_vector[sub_pop+8]->fitness, ga_vector[sub_pop+40]->pop, ga_vector[sub_pop+40]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+9]->pop, ga_vector[sub_pop+9]->fitness, ga_vector[sub_pop+41]->pop, ga_vector[sub_pop+41]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+10]->pop, ga_vector[sub_pop+10]->fitness, ga_vector[sub_pop+42]->pop, ga_vector[sub_pop+42]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+11]->pop, ga_vector[sub_pop+11]->fitness, ga_vector[sub_pop+43]->pop, ga_vector[sub_pop+43]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+12]->pop, ga_vector[sub_pop+12]->fitness, ga_vector[sub_pop+44]->pop, ga_vector[sub_pop+44]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+13]->pop, ga_vector[sub_pop+13]->fitness, ga_vector[sub_pop+45]->pop, ga_vector[sub_pop+45]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+14]->pop, ga_vector[sub_pop+14]->fitness, ga_vector[sub_pop+46]->pop, ga_vector[sub_pop+46]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+15]->pop, ga_vector[sub_pop+15]->fitness, ga_vector[sub_pop+47]->pop, ga_vector[sub_pop+47]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+16]->pop, ga_vector[sub_pop+16]->fitness, ga_vector[sub_pop+48]->pop, ga_vector[sub_pop+48]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+17]->pop, ga_vector[sub_pop+17]->fitness, ga_vector[sub_pop+49]->pop, ga_vector[sub_pop+49]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+18]->pop, ga_vector[sub_pop+18]->fitness, ga_vector[sub_pop+50]->pop, ga_vector[sub_pop+50]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+19]->pop, ga_vector[sub_pop+19]->fitness, ga_vector[sub_pop+51]->pop, ga_vector[sub_pop+51]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+20]->pop, ga_vector[sub_pop+20]->fitness, ga_vector[sub_pop+52]->pop, ga_vector[sub_pop+52]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+21]->pop, ga_vector[sub_pop+21]->fitness, ga_vector[sub_pop+53]->pop, ga_vector[sub_pop+53]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+22]->pop, ga_vector[sub_pop+22]->fitness, ga_vector[sub_pop+54]->pop, ga_vector[sub_pop+54]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+23]->pop, ga_vector[sub_pop+23]->fitness, ga_vector[sub_pop+55]->pop, ga_vector[sub_pop+55]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+24]->pop, ga_vector[sub_pop+24]->fitness, ga_vector[sub_pop+56]->pop, ga_vector[sub_pop+56]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+25]->pop, ga_vector[sub_pop+25]->fitness, ga_vector[sub_pop+57]->pop, ga_vector[sub_pop+57]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+26]->pop, ga_vector[sub_pop+26]->fitness, ga_vector[sub_pop+58]->pop, ga_vector[sub_pop+58]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+27]->pop, ga_vector[sub_pop+27]->fitness, ga_vector[sub_pop+59]->pop, ga_vector[sub_pop+59]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+28]->pop, ga_vector[sub_pop+28]->fitness, ga_vector[sub_pop+60]->pop, ga_vector[sub_pop+60]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+29]->pop, ga_vector[sub_pop+29]->fitness, ga_vector[sub_pop+61]->pop, ga_vector[sub_pop+61]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+30]->pop, ga_vector[sub_pop+30]->fitness, ga_vector[sub_pop+62]->pop, ga_vector[sub_pop+62]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+31]->pop, ga_vector[sub_pop+31]->fitness, ga_vector[sub_pop+63]->pop, ga_vector[sub_pop+63]->fitness, percentage);
            }

            // generations up to 210
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
                ga_vector[sub_pop]->evolve(30, false);
            }

            // exchange the best individuals
            for(int sub_pop = 0; sub_pop < no_sub_pops; sub_pop+=128) {
                exchange_percentage(ga_vector[sub_pop]->pop, ga_vector[sub_pop]->fitness, ga_vector[sub_pop+64]->pop, ga_vector[sub_pop+64]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+1]->pop, ga_vector[sub_pop+1]->fitness, ga_vector[sub_pop+65]->pop, ga_vector[sub_pop+65]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+2]->pop, ga_vector[sub_pop+2]->fitness, ga_vector[sub_pop+66]->pop, ga_vector[sub_pop+66]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+3]->pop, ga_vector[sub_pop+3]->fitness, ga_vector[sub_pop+67]->pop, ga_vector[sub_pop+67]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+4]->pop, ga_vector[sub_pop+4]->fitness, ga_vector[sub_pop+68]->pop, ga_vector[sub_pop+68]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+5]->pop, ga_vector[sub_pop+5]->fitness, ga_vector[sub_pop+69]->pop, ga_vector[sub_pop+69]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+6]->pop, ga_vector[sub_pop+6]->fitness, ga_vector[sub_pop+70]->pop, ga_vector[sub_pop+70]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+7]->pop, ga_vector[sub_pop+7]->fitness, ga_vector[sub_pop+71]->pop, ga_vector[sub_pop+71]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+8]->pop, ga_vector[sub_pop+8]->fitness, ga_vector[sub_pop+72]->pop, ga_vector[sub_pop+72]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+9]->pop, ga_vector[sub_pop+9]->fitness, ga_vector[sub_pop+73]->pop, ga_vector[sub_pop+73]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+10]->pop, ga_vector[sub_pop+10]->fitness, ga_vector[sub_pop+74]->pop, ga_vector[sub_pop+74]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+11]->pop, ga_vector[sub_pop+11]->fitness, ga_vector[sub_pop+75]->pop, ga_vector[sub_pop+75]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+12]->pop, ga_vector[sub_pop+12]->fitness, ga_vector[sub_pop+76]->pop, ga_vector[sub_pop+76]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+13]->pop, ga_vector[sub_pop+13]->fitness, ga_vector[sub_pop+77]->pop, ga_vector[sub_pop+77]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+14]->pop, ga_vector[sub_pop+14]->fitness, ga_vector[sub_pop+78]->pop, ga_vector[sub_pop+78]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+15]->pop, ga_vector[sub_pop+15]->fitness, ga_vector[sub_pop+79]->pop, ga_vector[sub_pop+79]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+16]->pop, ga_vector[sub_pop+16]->fitness, ga_vector[sub_pop+80]->pop, ga_vector[sub_pop+80]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+17]->pop, ga_vector[sub_pop+17]->fitness, ga_vector[sub_pop+81]->pop, ga_vector[sub_pop+81]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+18]->pop, ga_vector[sub_pop+18]->fitness, ga_vector[sub_pop+82]->pop, ga_vector[sub_pop+82]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+19]->pop, ga_vector[sub_pop+19]->fitness, ga_vector[sub_pop+83]->pop, ga_vector[sub_pop+83]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+20]->pop, ga_vector[sub_pop+20]->fitness, ga_vector[sub_pop+84]->pop, ga_vector[sub_pop+84]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+21]->pop, ga_vector[sub_pop+21]->fitness, ga_vector[sub_pop+85]->pop, ga_vector[sub_pop+85]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+22]->pop, ga_vector[sub_pop+22]->fitness, ga_vector[sub_pop+86]->pop, ga_vector[sub_pop+86]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+23]->pop, ga_vector[sub_pop+23]->fitness, ga_vector[sub_pop+87]->pop, ga_vector[sub_pop+87]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+24]->pop, ga_vector[sub_pop+24]->fitness, ga_vector[sub_pop+88]->pop, ga_vector[sub_pop+88]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+25]->pop, ga_vector[sub_pop+25]->fitness, ga_vector[sub_pop+89]->pop, ga_vector[sub_pop+89]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+26]->pop, ga_vector[sub_pop+26]->fitness, ga_vector[sub_pop+90]->pop, ga_vector[sub_pop+90]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+27]->pop, ga_vector[sub_pop+27]->fitness, ga_vector[sub_pop+91]->pop, ga_vector[sub_pop+91]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+28]->pop, ga_vector[sub_pop+28]->fitness, ga_vector[sub_pop+92]->pop, ga_vector[sub_pop+92]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+29]->pop, ga_vector[sub_pop+29]->fitness, ga_vector[sub_pop+93]->pop, ga_vector[sub_pop+93]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+30]->pop, ga_vector[sub_pop+30]->fitness, ga_vector[sub_pop+94]->pop, ga_vector[sub_pop+94]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+31]->pop, ga_vector[sub_pop+31]->fitness, ga_vector[sub_pop+95]->pop, ga_vector[sub_pop+95]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+32]->pop, ga_vector[sub_pop+32]->fitness, ga_vector[sub_pop+96]->pop, ga_vector[sub_pop+96]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+33]->pop, ga_vector[sub_pop+33]->fitness, ga_vector[sub_pop+97]->pop, ga_vector[sub_pop+97]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+34]->pop, ga_vector[sub_pop+34]->fitness, ga_vector[sub_pop+98]->pop, ga_vector[sub_pop+98]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+35]->pop, ga_vector[sub_pop+35]->fitness, ga_vector[sub_pop+99]->pop, ga_vector[sub_pop+99]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+36]->pop, ga_vector[sub_pop+36]->fitness, ga_vector[sub_pop+100]->pop, ga_vector[sub_pop+100]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+37]->pop, ga_vector[sub_pop+37]->fitness, ga_vector[sub_pop+101]->pop, ga_vector[sub_pop+101]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+38]->pop, ga_vector[sub_pop+38]->fitness, ga_vector[sub_pop+102]->pop, ga_vector[sub_pop+102]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+39]->pop, ga_vector[sub_pop+39]->fitness, ga_vector[sub_pop+103]->pop, ga_vector[sub_pop+103]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+40]->pop, ga_vector[sub_pop+40]->fitness, ga_vector[sub_pop+104]->pop, ga_vector[sub_pop+104]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+41]->pop, ga_vector[sub_pop+41]->fitness, ga_vector[sub_pop+105]->pop, ga_vector[sub_pop+105]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+42]->pop, ga_vector[sub_pop+42]->fitness, ga_vector[sub_pop+106]->pop, ga_vector[sub_pop+106]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+43]->pop, ga_vector[sub_pop+43]->fitness, ga_vector[sub_pop+107]->pop, ga_vector[sub_pop+107]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+44]->pop, ga_vector[sub_pop+44]->fitness, ga_vector[sub_pop+108]->pop, ga_vector[sub_pop+108]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+45]->pop, ga_vector[sub_pop+45]->fitness, ga_vector[sub_pop+109]->pop, ga_vector[sub_pop+109]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+46]->pop, ga_vector[sub_pop+46]->fitness, ga_vector[sub_pop+110]->pop, ga_vector[sub_pop+110]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+47]->pop, ga_vector[sub_pop+47]->fitness, ga_vector[sub_pop+111]->pop, ga_vector[sub_pop+111]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+48]->pop, ga_vector[sub_pop+48]->fitness, ga_vector[sub_pop+112]->pop, ga_vector[sub_pop+112]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+49]->pop, ga_vector[sub_pop+49]->fitness, ga_vector[sub_pop+113]->pop, ga_vector[sub_pop+113]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+50]->pop, ga_vector[sub_pop+50]->fitness, ga_vector[sub_pop+114]->pop, ga_vector[sub_pop+114]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+51]->pop, ga_vector[sub_pop+51]->fitness, ga_vector[sub_pop+115]->pop, ga_vector[sub_pop+115]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+52]->pop, ga_vector[sub_pop+52]->fitness, ga_vector[sub_pop+116]->pop, ga_vector[sub_pop+116]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+53]->pop, ga_vector[sub_pop+53]->fitness, ga_vector[sub_pop+117]->pop, ga_vector[sub_pop+117]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+54]->pop, ga_vector[sub_pop+54]->fitness, ga_vector[sub_pop+118]->pop, ga_vector[sub_pop+118]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+55]->pop, ga_vector[sub_pop+55]->fitness, ga_vector[sub_pop+119]->pop, ga_vector[sub_pop+119]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+56]->pop, ga_vector[sub_pop+56]->fitness, ga_vector[sub_pop+120]->pop, ga_vector[sub_pop+120]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+57]->pop, ga_vector[sub_pop+57]->fitness, ga_vector[sub_pop+121]->pop, ga_vector[sub_pop+121]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+58]->pop, ga_vector[sub_pop+58]->fitness, ga_vector[sub_pop+122]->pop, ga_vector[sub_pop+122]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+59]->pop, ga_vector[sub_pop+59]->fitness, ga_vector[sub_pop+123]->pop, ga_vector[sub_pop+123]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+60]->pop, ga_vector[sub_pop+60]->fitness, ga_vector[sub_pop+124]->pop, ga_vector[sub_pop+124]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+61]->pop, ga_vector[sub_pop+61]->fitness, ga_vector[sub_pop+125]->pop, ga_vector[sub_pop+125]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+62]->pop, ga_vector[sub_pop+62]->fitness, ga_vector[sub_pop+126]->pop, ga_vector[sub_pop+126]->fitness, percentage);
                exchange_percentage(ga_vector[sub_pop+63]->pop, ga_vector[sub_pop+63]->fitness, ga_vector[sub_pop+127]->pop, ga_vector[sub_pop+127]->fitness, percentage);
            }
        }

        for (int sub_pop = 0; sub_pop < no_sub_pops; sub_pop++) {
            ga_vector[sub_pop]->~Ga();
        }
    }


# endif

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    printf("Time taken by function: %lld seconds", duration.count());
    return 0;
}