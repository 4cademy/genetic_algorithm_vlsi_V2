//
// Created by Marcel Beyer on 11.09.2023.
//

#ifndef GENETIC_ALGORITHM_VLSI_V2_GA_H
#define GENETIC_ALGORITHM_VLSI_V2_GA_H


class Ga {
public:
    Ga()= default;;
    Ga(unsigned dim, unsigned pop_size, double min_gene, double max_gene);
    void load_data_f1();
    void evolve(int generations);

    void clean() const;

    double** pop{};
    double* fitness{};
    double min_fitness{};
    double max_fitness{};


private:
    unsigned dim{};
    unsigned pop_size{};
    double min_gene{};
    double max_gene{};
    double* opt1{};
    double** mating_list{};  // contains the numbers of the individuals that are selected for mating (index0 mates with index1, index2 with index3, ...)

    double function1(const double* individual) const;
    void compute_fitness();
    void selection_roulette();
    void crossover_uniform();
};


#endif //GENETIC_ALGORITHM_VLSI_V2_GA_H
