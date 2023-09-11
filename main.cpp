//
// Created by Marcel Beyer on 22.08.2023.
//
#include <random>
#include "population.h"
#include "ga_old.h"
#include <chrono>
using namespace std::chrono;
#include <limits>

unsigned dim = 1000;
unsigned pop_size = 5000;

std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<double> uniformRealDistribution(0.0, 1.0);

int main() {
    auto start = high_resolution_clock::now();
    Population population0(dim,pop_size, -100, 100);
    Population population1(dim,pop_size, -100, 100);

    GaOld ga0(population0, dim);
    GaOld ga1(population1, dim);

    for(int i = 0; i < 10; i++) {
        printf("Generation %d\n", i);
        printf("%e\n", population0.fitness.min_fitness);
        printf("%e\n", population1.fitness.min_fitness);
        ga0.evolve(100);  // -1 for unlimited generations until convergence is reached
        ga1.evolve(100);  // -1 for unlimited generations until convergence is reached
        auto** min_list = new double*[pop_size];
        for (unsigned j = 0; j < pop_size; j++) {
            min_list[j] = new double[dim];
        }
        for(int j = 0; j < pop_size; j++) {
            int min0 = 0;
            int min1 = 0;
            for(int k = 0; k < pop_size; k++) {
                if(population0.fitness.fitness_values[k] < population0.fitness.fitness_values[min0]) {
                    min0 = k;
                }
                if(population1.fitness.fitness_values[k] < population1.fitness.fitness_values[min1]) {
                    min1 = k;
                }
                if(population0.fitness.fitness_values[min0] < population1.fitness.fitness_values[min1]) {
                    for(int l = 0; l < dim; l++) {
                        min_list[j][l] = population0.pop[min0][l];
                        population0.fitness.fitness_values[min0] = std::numeric_limits<double>::max();
                    }
                } else {
                    for(int l = 0; l < dim; l++) {
                        min_list[j][l] = population1.pop[min1][l];
                        population1.fitness.fitness_values[min1] = std::numeric_limits<double>::max();
                    }
                }
            }
        }
        for(int j = 0; j < pop_size; j++) {
            for(int k = 0; k < dim; k++) {
                population0.pop[j][k] = min_list[j][k];
                population1.pop[j][k] = min_list[j][k];
            }
        }
    }

    // free allocated memory
    population0.clean();
    population1.clean();
    ga0.clean();
    ga1.clean();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    printf("Time taken by function: %lld microseconds", duration.count());
    return 0;
}