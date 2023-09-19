//
// Created by Marcel Beyer on 11.09.2023.
//

#include <cstdio>
#include <random>
#include "ga.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <omp.h>
#include <cmath>
#include <thread>


std::random_device rd;
const std::random_device::result_type seed = rd();

// loads the data for the optimal solution for function 1 from the datafile into the opt1 vector
void Ga::load_data_f1(){
    opt1 = new double[dim];
    int i = 0;
    std::string path = "../cdatafiles/F1-xopt.txt";
    std::ifstream file(path);
    if (file.is_open()) {
        std::string line;
        while (getline(file, line)) {
            opt1[i++] = stod(line);
        }
        file.close();
    }
    else {
        std::cout << "Cannot open the datafile '" << path << "'" << std::endl;
    }
    printf("%f\n", opt1[0]);
}

// fitness function 1
double Ga::function1(const double* individual) const {
    double result = 0;
    auto* z = new double[dim];

    int sign;
    double hat;
    double c1;
    double c2;

    # pragma omp parallel for default(none) shared(individual, z, dim, opt1) private(sign, hat, c1, c2) reduction(+:result)
    for(unsigned i = 0; i < dim; i++) {
        z[i] = individual[i] - opt1[i];
        // Transformation
        if (z[i] > 0) {
            sign = 1;
            hat = log(fabs(z[i]));
            c1 = 10;
            c2 = 7.9;
            z[i] = sign*exp(hat+0.049*(sin(c1*hat)+sin(c2*hat)));
        } else if (z[i] == 0) {
            z[i] = 0;
        } else {
            sign = -1;
            hat = log(fabs(z[i]));
            c1 = 5.5;
            c2 = 3.1;
            z[i] = sign*exp(hat+0.049*(sin(c1*hat)+sin(c2*hat)));
        }
        result += pow(1.0e6,  i/((double)(dim - 1)) ) * z[i] * z[i];
    }

    delete[] z;
    return result;
}

// computes the fitness of each individual in the population
void Ga::compute_fitness(){
    fitness[0] = function1(pop[0]);
    min_fitness = function1(pop[0]);
    max_fitness = function1(pop[0]);
    #pragma omp parallel for default(none) shared(pop, fitness, pop_size) reduction(min:min_fitness) reduction(max:max_fitness)
    for (unsigned i=1; i<pop_size; i++) {
        fitness[i] = function1(pop[i]);
        if(fitness[i] < min_fitness) {
            min_fitness = fitness[i];
        } else if (fitness[i] > max_fitness) {
            max_fitness = fitness[i];
        }
    }

    // checks if the current min_fitness is better than the best_fitness
    if (min_fitness < best_fitness) {
        best_fitness = min_fitness;
    }
    // computes the convergence
    convergence = 1-min_fitness/max_fitness;
}

// constructor for the class Ga
Ga::Ga(unsigned int dim, unsigned int pop_size, double min_gene, double max_gene) {
    this->dim = dim;
    this->pop_size = pop_size;
    this->min_gene = min_gene;
    this->max_gene = max_gene;

    // load optimal vector for benchmark function 1
    load_data_f1();

    // initialize population
    pop = new double*[pop_size];
    for (unsigned i = 0; i < pop_size; i++) {
        pop[i] = new double[dim];
    }

    #pragma omp parallel default(none) shared(seed, pop_size, dim, pop, min_gene, max_gene)
    {
        std::uniform_real_distribution<double> dist(Ga::min_gene, Ga::max_gene);
        std::mt19937_64 gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + clock());
        #pragma omp for collapse(2)
        for (unsigned i = 0; i < pop_size; i++) {
            for (unsigned j = 0; j < dim; j++) {
                pop[i][j] = dist(gen);
            }
        }
    }

    // initialize fitness vector
    fitness = new double[pop_size];
    compute_fitness();
    best_fitness = min_fitness;

    min_fitness_vector = new double[pop_size];

    // initialize mating list
    mating_list = new double*[2*pop_size];
    for (unsigned i = 0; i < 2*pop_size; i++) {
        mating_list[i] = new double[dim];
    }
}

// cleans up the allocated memory
void Ga::clean() const {
    for (unsigned i = 0; i < pop_size; i++) {
        delete[] pop[i];
    }
    delete[] pop;

    for (unsigned i = 0; i < pop_size; i++) {
        delete[] mating_list[i];
    }
    delete[] mating_list;

    delete[] opt1;
    delete[] fitness;
}

void Ga::selection_roulette() {
    // create list of offsets for roulette selection
    auto* offset = new double[pop_size];
    auto* roulette_fitness = new double[pop_size]; // an all positive fitness vector where the best individual has the highest fitness
    double total_fitness = 0;

    // shift fitness to all positive values and invert it so the minimal value has the highest fitness
    #pragma omp parallel for default(none) shared(pop_size, roulette_fitness, max_fitness, fitness) reduction(+:total_fitness)
    for (unsigned i=0; i<pop_size; i++) {
        // shift fitness to all positive values and invert it so the minimal value has the highest fitness
        roulette_fitness[i] = max_fitness - fitness[i];
        // calculate total fitness
        total_fitness += roulette_fitness[i];
    }

    // ToDo: merge with for-loop for roulette_fitness
    offset[0] = roulette_fitness[0] / total_fitness;
    // calculate offset for roulette selection
    for (unsigned i=1; i<pop_size; i++) {
        offset[i] = offset[i-1] + (roulette_fitness[i] / total_fitness);
    }

    // create global random number generator for use outside of OpenMP parallel regions
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::mt19937_64 gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + clock());

    // do roulette selection
    for (unsigned i=0; i < 2*pop_size; i++) {
        double roulette_random = dist(gen);
        for (unsigned j=0; j<pop_size; j++) {
            if (roulette_random < offset[j]) {
                for (unsigned k=0; k < dim; k++) {
                    mating_list[i][k] = pop[j][k];          // ToDo prevent mating with itself
                }
                break;
            }
        }
    }

    // free memory
    delete[] offset;
    delete[] roulette_fitness;
}

void Ga::crossover_uniform() {
    #pragma omp parallel default(none) shared(seed, pop_size, dim, pop, mating_list)
    {
        std::uniform_real_distribution<double> dist(0, 1);
        std::mt19937_64 gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + clock());
        #pragma for collapse(2)
        for (unsigned i = 0; i < pop_size; i++) {
            for (unsigned j = 0; j < dim; j++) {
                if (dist(gen) < 0.5) {       // choose gene of parent A
                    pop[i][j] = mating_list[2 * i][j];
                } else {                    // choose gene of parent B
                    pop[i][j] = mating_list[2 * i + 1][j];
                }
            }
        }
    }
}

void Ga::evolve(int generations, bool break_on_convergence) {
    int i;
    unsigned convergence_counter = 0;
    double convergence_threshold = 0.1;
    for (i=0; i<generations; i++) {
        compute_fitness();
        Ga::min_fitness_vector[i] = min_fitness;

        selection_roulette();
        crossover_uniform();
        // ToDo mutation
        // printf("%e \n", min_fitness);
        if(break_on_convergence){
            if (convergence < convergence_threshold) {
                convergence_counter++;
                if (convergence_counter >= 10){
                    break;
                }
            } else {
                convergence_counter = 0;
            }
        }
    }

    // write results to file
    const auto p1 = std::chrono::system_clock::now();
    long long timestamp = std::chrono::duration_cast<std::chrono::seconds>(p1.time_since_epoch()).count();
    std::ofstream results_file;
    std::string filename = "results_" + std::to_string(timestamp) + ".csv";
    results_file.open(filename);
    for (unsigned j = 0; j < i; j++) {
        results_file << std::to_string(Ga::min_fitness_vector[j]) << "\n";
    }

    printf("Generations: %d \n", i);
    printf("Best fitness: %e \n", Ga::best_fitness);
}


