//
// Created by Marcel Beyer on 11.09.2023.
//

#include <vector>

#ifndef GENETIC_ALGORITHM_VLSI_V2_GA_H
#define GENETIC_ALGORITHM_VLSI_V2_GA_H


class Ga {
public:
    Ga()= delete;
    Ga(unsigned dim, unsigned pop_size, float min_gene, float max_gene);
    Ga(const Ga& ga) = delete;
    Ga(Ga&& ga) = delete;
    ~Ga();
    void load_data_f1();
    void evolve(int generations, bool break_on_convergence);

    unsigned dim{};
    unsigned pop_size{};
    float min_gene{};  // min value for each gene
    float max_gene{};  // max value for each gene
    float** pop{};
    float* fitness{};
    float min_fitness{};
    unsigned min_fitness_index{};
    float max_fitness{};
    unsigned max_fitness_index{};
    float best_fitness{};
    float convergence{};
    float convergence_counter{};
    float mutation_rate{};
    float mutation_deviation{};


private:
    float* opt1{};
    float** mating_list{};  // contains the numbers of the individuals that are selected for mating (index0 mates with index1, index2 with index3, ...)
    // contains the minimum fitness of each generation
    long long start_time{}; // start time of the algorithm

    void selection_roulette();
    void crossover_uniform();
    void mutation_normal_dist();

protected:
    std::vector<float> min_fitness_vector;

    void reset_population();

    float function1(const float* individual) const;

    void compute_fitness();
};


#endif //GENETIC_ALGORITHM_VLSI_V2_GA_H
