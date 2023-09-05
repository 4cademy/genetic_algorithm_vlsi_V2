//
// Created by Marcel Beyer on 05.09.2023.
//

#ifndef GENETIC_ALGORITHM_VLSI_V2_POPULATION_H
#define GENETIC_ALGORITHM_VLSI_V2_POPULATION_H


class Population {
public:
    explicit Population(unsigned pop_size, double min, double max);
    void debug_print() const;
    void clean() const;

private:
    unsigned pop_size;
};


#endif //GENETIC_ALGORITHM_VLSI_V2_POPULATION_H
