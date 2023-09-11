//
// Created by Marcel Beyer on 04.09.2023.
//

#ifndef GENETIC_ALGORITHM_VLSI_V2_FITNESS_H
#define GENETIC_ALGORITHM_VLSI_V2_FITNESS_H


class Fitness {
public:
    Fitness()= default;;
    Fitness(double** pop, unsigned pop_size);
    static void debug_print();
    static double function1(const double *individual);
    static void clean();
    double* fitness_values{};
    double min_fitness{};
    double max_fitness{};
    unsigned pop_size;
    void compute_fitness(double **pop);
};


#endif //GENETIC_ALGORITHM_VLSI_V2_FITNESS_H