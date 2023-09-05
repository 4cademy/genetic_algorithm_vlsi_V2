//
// Created by Marcel Beyer on 04.09.2023.
//

#ifndef GENETIC_ALGORITHM_VLSI_V2_FITNESS_H
#define GENETIC_ALGORITHM_VLSI_V2_FITNESS_H


class Fitness {
public:
    explicit Fitness(double** pop, unsigned pop_size);
    static void debug_print();
    double get_fitness(unsigned index);
    double function1(const double *individual);
    static void clean();
private:
    double* fitness_values;
};


#endif //GENETIC_ALGORITHM_VLSI_V2_FITNESS_H