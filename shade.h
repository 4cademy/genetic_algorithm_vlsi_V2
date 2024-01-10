//
// Created by Marcel Beyer on 10.01.2024.
//

#ifndef GENETIC_ALGORITHM_VLSI_V2_SHADE_H
#define GENETIC_ALGORITHM_VLSI_V2_SHADE_H

#include "ga.h"

class ShadeGa : Ga {
public:
    ShadeGa()= delete;
    ShadeGa(unsigned dim, unsigned pop_size, float min_gene, float max_gene);
    ShadeGa(const ShadeGa& ga) = delete;
    ShadeGa(ShadeGa&& ga) = delete;
    ~ShadeGa();
    void evolve_shade(int generations);

    float prev_min_fitness{};

private:
    unsigned history_size{};
    unsigned history_replacement_index{};
    float* cr_array{};
    float* f_array{};
    float* s_cr{};
    float* s_f{};
    float* s_delta_fitness{};
    unsigned s_index{};
    float** trial_pop{};
    float* trial_fitness{};
    float* trial_cr{};
    float* trial_f{};
    float rate_of_improvement{};
    void generate_trial_vector(float* trial_vector, unsigned parent_vector_index, float &return_cr, float &return_f);
    void shade();
};

#endif //GENETIC_ALGORITHM_VLSI_V2_SHADE_H
