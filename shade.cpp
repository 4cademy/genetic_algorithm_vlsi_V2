//
// Created by marce on 10.01.2024.
//

#include <cstdio>
#include <cmath>
#include <random>
#include <thread>
#include "shade.h"

std::random_device shade_rd;
const std::random_device::result_type seed = shade_rd();

ShadeGa::ShadeGa(unsigned int dim, unsigned int pop_size, float min_gene, float max_gene) : Ga(dim, pop_size, min_gene, max_gene) {
    // initialize shade
    history_size = 100;
    history_replacement_index = 0;
    cr_array = new float[history_size];
    f_array = new float[history_size];
    s_cr = new float[pop_size];
    s_f = new float[pop_size];
    s_delta_fitness = new float[pop_size];
    s_index = 0;
    for (unsigned i = 0; i < history_size; i++) {
        cr_array[i] = 0.5;
        f_array[i] = 0.5;
    }
    trial_pop = new float*[pop_size];
    for (int i = 0; i < pop_size; i++) {
        trial_pop[i] = new float[dim];
    }
    trial_fitness = new float[pop_size];
    trial_cr = new float[pop_size];
    trial_f = new float[pop_size];
}

ShadeGa::~ShadeGa() {
    delete[] cr_array;
    delete[] f_array;
    delete[] s_cr;
    delete[] s_f;
    delete[] s_delta_fitness;
    for (int i = 0; i < pop_size; i++) {
        delete[] trial_pop[i];
    }
    delete[] trial_pop;
    delete[] trial_fitness;
    delete[] trial_cr;
    delete[] trial_f;
}

float weighted_arithmetic_mean(const float* values, const float* weights, unsigned max_index){
    float total_weight = 0;
    for(int i = 0; i <= max_index; i++){
        total_weight += weights[i];
    }

    float mean = 0;
    for(int i = 0; i <= max_index; i++){
        mean += (values[i] * weights[i]/total_weight);
    }
    return mean;
}

float weighted_squared_mean(const float* values, const float* weights, unsigned max_index){
    float total_weight = 0;
    for(int i = 0; i <= max_index; i++){
        total_weight += weights[i];
    }

    float mean = 0;
    for(int i = 0; i <= max_index; i++){
        mean += (values[i] * values[i] * weights[i]/total_weight);
    }
    return mean;
}

float weighted_lehmer_mean(const float* values, const float* weights, unsigned max_index){
    float mean_of_squares = weighted_squared_mean(values, weights, max_index);
    float arithmetic_mean = weighted_arithmetic_mean(values, weights, max_index);
    float mean = mean_of_squares/arithmetic_mean;
    return mean;
}

void ShadeGa::generate_trial_vector(float* trial_vector, unsigned parent_vector_index, float &return_cr, float &return_f) {
    const auto p1 = std::chrono::system_clock::now();
    long long timestamp = std::chrono::duration_cast<std::chrono::microseconds>(p1.time_since_epoch()).count();

    std::uniform_int_distribution<unsigned> uni_int_dist(0, history_size-1);
    std::mt19937_64 uni_int_gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + timestamp);

    unsigned r1 = uni_int_dist(uni_int_gen);
    float cr = cr_array[r1];
    float f = f_array[r1];

    std::normal_distribution<float> normal_float_dist(cr, 0.1);
    std::mt19937_64 normal_float_gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + timestamp + 1);

    std::cauchy_distribution<float> cauchy_dist(f, 0.1);
    std::mt19937_64 cauchy_gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + timestamp + 2);

    std::uniform_real_distribution<float> uni_dist(0, 1);
    std::mt19937_64 uni_gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + timestamp + 3);

    cr = normal_float_dist(normal_float_gen);
    if (cr < 0) {
        cr = 0;
    } else if (cr > 1) {
        cr = 1;
    }

    f = cauchy_dist(cauchy_gen);
    while (f < 0) {
        f = cauchy_dist(cauchy_gen);
    }
    if (f > 1) {
        f = 1;
    }

    // ToDo: implement current-to-pbest/1/bin
    uni_int_dist = std::uniform_int_distribution<unsigned>(0, pop_size-1);
    unsigned j = uni_int_dist(uni_int_gen);
    unsigned donor1_index = uni_int_dist(uni_int_gen);
    unsigned donor2_index = uni_int_dist(uni_int_gen);

    while (donor1_index == parent_vector_index) {
        donor1_index = uni_int_dist(uni_int_gen);
    }

    while ((donor2_index == donor1_index) || (donor2_index == parent_vector_index)) {
        donor2_index = uni_int_dist(uni_int_gen);
    }

    for (unsigned i = 0; i < dim; i++) {
        float old_gene = pop[parent_vector_index][i];
        float diff = f * ((pop[Ga::min_fitness_index][i] - old_gene) + (pop[donor1_index][i] - pop[donor2_index][i]));
        float new_gene = old_gene + diff;
        if (new_gene < min_gene) {
            new_gene = (min_gene + new_gene)/2;
        } else if (new_gene > max_gene) {
            new_gene = (max_gene + new_gene)/2;
        }

        if ((j == i) || (uni_dist(uni_gen) < cr)) {
            trial_vector[i] = new_gene;
        } else {
            trial_vector[i] = old_gene;
        }
    }

    return_cr = cr;
    return_f = f;
}

void ShadeGa::shade() {
    ShadeGa::prev_min_fitness = Ga::min_fitness;

    for (int i = 0; i < pop_size; i++) {
        generate_trial_vector(trial_pop[i], i, trial_cr[i], trial_f[i]);
        trial_fitness[i] = Ga::function1(trial_pop[i]);
    }

    for (int i = 0; i < pop_size; i++) {
        if (trial_fitness[i] < fitness[i]) {
            s_cr[s_index] = trial_cr[i];
            s_f[s_index] = trial_f[i];
            s_delta_fitness[s_index] = fabsf(fitness[i]-trial_fitness[i]);
            s_index++;

            for (int j = 0; j < dim; j++) {
                pop[i][j] = trial_pop[i][j];
            }
            fitness[i] = trial_fitness[i];
            if (trial_fitness[i] < min_fitness) {
                min_fitness = trial_fitness[i];
                min_fitness_index = i;
            }
            if (trial_fitness[i] < best_fitness) {
                best_fitness = trial_fitness[i];
            }
        }
    }

    float mean_wa_cr = 0;
    float mean_wl_f = 0;
    if (s_index > 0){
        mean_wa_cr = weighted_arithmetic_mean(s_cr, s_delta_fitness, s_index-1);       // weighted arithmetic mean
        mean_wl_f = weighted_lehmer_mean(s_cr, s_delta_fitness, s_index-1);        // weighted Lehmer mean
        cr_array[history_replacement_index] = mean_wa_cr;
        f_array[history_replacement_index] = mean_wl_f;
        history_replacement_index++;
        if (history_replacement_index >= history_size) {
            history_replacement_index = 0;
        }

        s_index = 0;
    }

    rate_of_improvement = (ShadeGa::prev_min_fitness - Ga::min_fitness) / ShadeGa::prev_min_fitness;
    printf("Rate of improvement: %f\n", rate_of_improvement);
    if (rate_of_improvement < 0) {
        Ga::convergence_counter++;
    } else {
        Ga::convergence_counter = 0;
    }
}

void ShadeGa::evolve_shade(int generations) {
    int i;
    for (i=0; i<generations; i++) {
        shade();

        Ga::min_fitness_vector.push_back(min_fitness);

        printf("%i: %e\n", i, min_fitness);

        if (Ga::convergence_counter >= 100) {
            reset_population();
            compute_fitness();

            Ga::convergence_counter = 0;
        }
    }

    printf("Generations: %d \n", i);
    printf("Best fitness: %e \n", Ga::best_fitness);
}
