//
// Created by Marcel Beyer on 04.09.2023.
//

#ifndef GENETIC_ALGORITHM_VLSI_V2_FITNESS_H
#define GENETIC_ALGORITHM_VLSI_V2_FITNESS_H
class Fitness {
public:
    explicit Fitness();
    static void debug_print();
    double function1(const double *x);
};
#endif //GENETIC_ALGORITHM_VLSI_V2_FITNESS_H