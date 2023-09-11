//
// Created by Marcel Beyer on 11.09.2023.
//

#include <random>
#include "ga.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <omp.h>
#include <cmath>

std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<double> uniformRealDistribution(0.0, 1.0);

// creates random number between min and max
double rng(double min, double max) {
    return uniformRealDistribution(gen) * (max-min) + min;
}

// loads the data for the optimal solution for function 1 from the datafile into the opt1 vector
void Ga::load_data_f1(){
    opt1 = new double[dim];
    int i = 0;
    std::string path = "../cdatafiles/F1-xopt.txt";
    std::ifstream file(path);
    printf("Debug: 4.2\n");
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

    // # pragma omp parallel for default(none) shared(individual, z, dim, opt1) private(sign, hat, c1, c2) reduction(+:result)
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
    for (unsigned i=0; i<pop_size; i++) {
        fitness[i] = function1(pop[i]);
    }
}

// constructor for the class Ga
Ga::Ga(unsigned int dim, unsigned int pop_size, double min_gene, double max_gene) {
    this->dim = dim;
    this->pop_size = pop_size;
    this->min_gene = min_gene;
    this->max_gene = max_gene;

    // initialize population
    pop = new double*[pop_size];
    for (unsigned i = 0; i < pop_size; i++) {
        pop[i] = new double[dim];
    }
    for (unsigned i=0; i<pop_size; i++) {
        for (unsigned j=0; j<dim; j++) {
            pop[i][j] = rng(min_gene, max_gene);
        }
    }

    // load optimal vector for benchmark function 1
    load_data_f1();

    // initialize fitness vector
    fitness = new double[pop_size];
    compute_fitness();

}

// cleans up the allocated memory
void Ga::clean() const {
    for (unsigned i = 0; i < pop_size; i++) {
        delete[] pop[i];
    }
    delete[] pop;

    delete[] opt1;
    delete[] fitness;
}


