//
// Created by Marcel Beyer on 04.09.2023.
//

#ifndef GENETIC_ALGORITHM_VLSI_V2_FITNESS_H
#define GENETIC_ALGORITHM_VLSI_V2_FITNESS_H
class Fitness {
public:
    explicit Fitness(unsigned dimension);
private:
    static unsigned dim;
    static double* opt1;
};
#endif //GENETIC_ALGORITHM_VLSI_V2_FITNESS_H