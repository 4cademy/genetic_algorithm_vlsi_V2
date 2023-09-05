//
// Created by Marcel Beyer on 22.08.2023.
//
#include <random>
#include "population.h"

unsigned dim = 1000;

std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<double> uniformRealDistribution(0.0, 1.0);

int main() {
    Population population(dim,2, -100, 100);
    population.debug_print();

    // free allocated memory
    population.clean();

    return 0;
}