//
// Created by Marcel Beyer on 11.09.2023.
//

#include <cstdio>
#include <random>
#include "ga.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <omp.h>
#include <cmath>
#include <thread>
#include <vector>


std::random_device rd;
const std::random_device::result_type seed = rd();

// loads the data for the optimal solution for function 1 from the datafile into the opt1 vector
void Ga::load_data_f1(){
    opt1 = new float[dim];
    int i = 0;
    std::string path = "../cdatafiles/F1-xopt.txt";
    std::ifstream file(path);
    if (file.is_open()) {
        std::string line;
        while (getline(file, line)) {
            opt1[i++] = stof(line);
        }
        file.close();
    }
    else {
        std::cout << "Cannot open the datafile '" << path << "'" << std::endl;
    }
}

// fitness function 1
float Ga::function1(const float* individual) const {
    float result = 0;
    auto* z = new float[dim];

    int sign;
    float hat;
    float c1;
    float c2;

    # pragma omp parallel for default(none) shared(individual, z, dim, opt1) private(sign, hat, c1, c2) reduction(+:result)
    for(unsigned i = 0; i < dim; i++) {
        z[i] = individual[i] - opt1[i];
        // Transformation
        if (z[i] > 0) {
            sign = 1;
            hat = logf(fabsf(z[i]));
            c1 = 10;
            c2 = 7.9;
            z[i] = (float)sign*expf(hat+(float)0.049*(sinf(c1*hat)+sinf(c2*hat)));
        } else if (z[i] == 0) {
            z[i] = 0;
        } else {
            sign = -1;
            hat = logf(fabsf(z[i]));
            c1 = 5.5;
            c2 = 3.1;
            z[i] = (float)sign*expf(hat+(float)0.049*(sinf(c1*hat)+sinf(c2*hat)));
        }
        result += powf(1.0e6,  (float)i/(float)((dim - 1)) ) * z[i] * z[i];
    }

    delete[] z;
    return result;
}

// computes the fitness of each individual in the population
void Ga::compute_fitness(){
    fitness[0] = function1(pop[0]);
    min_fitness = function1(pop[0]);
    min_fitness_index = 0;
    max_fitness = function1(pop[0]);
    max_fitness_index = 0;
    #pragma omp parallel for default(none) shared(pop, fitness, pop_size) reduction(min:min_fitness) reduction(max:max_fitness)
    for (unsigned i=1; i<pop_size; i++) {
        fitness[i] = function1(pop[i]);
        if(fitness[i] < min_fitness) {
            min_fitness = fitness[i];
            min_fitness_index = i;
        } else if (fitness[i] > max_fitness) {
            max_fitness = fitness[i];
            max_fitness_index = i;
        }
    }

    // checks if the current min_fitness is better than the best_fitness
    if (min_fitness < best_fitness) {
        best_fitness = min_fitness;
    }
    // computes the convergence
    convergence = 1-min_fitness/max_fitness;
}

// constructor for the class Ga
Ga::Ga(unsigned int dim, unsigned int pop_size, float min_gene, float max_gene) {
    this->dim = dim;
    this->pop_size = pop_size;
    this->min_gene = min_gene;
    this->max_gene = max_gene;

    // load optimal vector for benchmark function 1
    load_data_f1();

    // initialize population
    pop = new float*[pop_size];
    for (unsigned i = 0; i < pop_size; i++) {
        pop[i] = new float[dim];
    }

    #pragma omp parallel default(none) shared(seed, pop_size, dim, pop, min_gene, max_gene)
    {
        const auto p1 = std::chrono::system_clock::now();
        long long timestamp = std::chrono::duration_cast<std::chrono::microseconds>(p1.time_since_epoch()).count();
        std::uniform_real_distribution<float> dist(Ga::min_gene, Ga::max_gene);
        std::mt19937_64 gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + timestamp);
        #pragma omp for collapse(2)
        for (unsigned i = 0; i < pop_size; i++) {
            for (unsigned j = 0; j < dim; j++) {
                pop[i][j] = dist(gen);
            }
        }
    }

    // initialize fitness vector
    fitness = new float[pop_size];
    compute_fitness();
    best_fitness = min_fitness;

    // initialize mating list
    mating_list = new float*[2*pop_size];
    for (unsigned i = 0; i < 2*pop_size; i++) {
        mating_list[i] = new float[dim];
    }

    // initialize mutation rate and deviation
    mutation_rate = 0.01;
    mutation_deviation = (max_gene - min_gene) / 6;

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

    // get start time
    const auto p1 = std::chrono::system_clock::now();
    Ga::start_time = std::chrono::duration_cast<std::chrono::microseconds>(p1.time_since_epoch()).count();

    //
    printf("Ga created: %p\n", this);
}

// destructor for the class Ga
Ga::~Ga() {
    // write results to file
    std::ofstream results_file;
    std::string filename = "results_" + std::to_string(Ga::start_time) + ".csv";
    results_file.open(R"(C:\Users\marce\OneDrive\Dokumente\ga_results\)" + filename, std::ios_base::app);
    for (float val : Ga::min_fitness_vector) {
        results_file << std::to_string(val) << "\n";
    }
    results_file.close();

    // free memory
    for (unsigned i = 0; i < pop_size; i++) {
        delete[] pop[i];
    }
    delete[] pop;
    for (unsigned i = 0; i < pop_size; i++) {
        delete[] mating_list[i];
    }
    delete[] mating_list;
    delete[] opt1;
    delete[] fitness;

    // free memory for shade
    for (int i = 0; i < pop_size; i++) {
        delete[] trial_pop[i];
    }
    delete[] trial_pop;
    delete[] trial_cr;
    delete[] trial_f;
    delete[] trial_fitness;
    printf("Ga destroyed: %p\n", this);
}


void Ga::selection_roulette() {
    // create list of offsets for roulette selection
    auto* offset = new float[pop_size];
    auto* roulette_fitness = new float[pop_size]; // an all positive fitness vector where the best individual has the highest fitness
    float total_fitness = 0;

    // shift fitness to all positive values and invert it so the minimal value has the highest fitness
    #pragma omp parallel for default(none) shared(pop_size, roulette_fitness, max_fitness, fitness) reduction(+:total_fitness)
    for (unsigned i=0; i<pop_size; i++) {
        // shift fitness to all positive values and invert it so the minimal value has the highest fitness
        roulette_fitness[i] = max_fitness - fitness[i];
        // calculate total fitness
        total_fitness += roulette_fitness[i];
    }

    // ToDo: merge with for-loop for roulette_fitness
    offset[0] = roulette_fitness[0] / total_fitness;
    // calculate offset for roulette selection
    for (unsigned i=1; i<pop_size; i++) {
        offset[i] = offset[i-1] + (roulette_fitness[i] / total_fitness);
    }

    // create global random number generator for use outside of OpenMP parallel regions
    const auto p1 = std::chrono::system_clock::now();
    long long timestamp = std::chrono::duration_cast<std::chrono::microseconds>(p1.time_since_epoch()).count();
    std::uniform_real_distribution<float> dist(0.0, 1.0);
    std::mt19937_64 gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + timestamp);
    // ToDo: parallelize
    // do roulette selection
    for (unsigned i=0; i < 2*pop_size; i++) {
        float roulette_random = dist(gen);
        for (unsigned j=0; j<pop_size; j++) {
            if (roulette_random < offset[j]) {
                for (unsigned k=0; k < dim; k++) {
                    mating_list[i][k] = pop[j][k];          // ToDo prevent mating with itself
                }
                break;
            }
        }
    }

    // free memory
    delete[] offset;
    delete[] roulette_fitness;
}

void Ga::crossover_uniform() {
    #pragma omp parallel default(none) shared(seed, pop_size, dim, pop, mating_list)
    {
        const auto p1 = std::chrono::system_clock::now();
        long long timestamp = std::chrono::duration_cast<std::chrono::microseconds>(p1.time_since_epoch()).count();
        std::uniform_real_distribution<float> dist(0, 1);
        std::mt19937_64 gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + timestamp);
        #pragma for collapse(2)
        for (unsigned i = 0; i < pop_size; i++) {
            for (unsigned j = 0; j < dim; j++) {
                if (dist(gen) < 0.5) {       // choose gene of parent A
                    pop[i][j] = mating_list[2 * i][j];
                } else {                    // choose gene of parent B
                    pop[i][j] = mating_list[2 * i + 1][j];
                }
            }
        }
    }
}

void Ga::mutation_normal_dist() {
    const auto p1 = std::chrono::system_clock::now();
    long long timestamp = std::chrono::duration_cast<std::chrono::microseconds>(p1.time_since_epoch()).count();

    std::normal_distribution<float> normal_dist(0, mutation_deviation);    // standard deviation is 1/6 of the range -> 99.8% of the values are in the range
    std::mt19937_64 normal_gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + timestamp);

    std::uniform_real_distribution<float> uni_dist(0, 1);
    std::mt19937_64 uni_gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + timestamp + 1);

    for(unsigned i = 0; i < pop_size; i++) {
        for(unsigned j = 0; j < dim; j++) {
            if(uni_dist(uni_gen) < mutation_rate) {
                float mutated_gene = pop[i][j] + normal_dist(normal_gen);
                while(mutated_gene < min_gene || mutated_gene > max_gene) {
                    mutated_gene = pop[i][j] + normal_dist(normal_gen);
                }
                pop[i][j] = mutated_gene;
            }
        }
    }

}

void Ga::generate_trial_vector(float* trial_vector, unsigned parent_vector_index, float &return_cr, float &return_f) {
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

void Ga::reset_population() {
    const auto p1 = std::chrono::system_clock::now();
    long long timestamp = std::chrono::duration_cast<std::chrono::microseconds>(p1.time_since_epoch()).count();
    std::uniform_real_distribution<float> dist(Ga::min_gene, Ga::max_gene);
    std::mt19937_64 gen(seed + std::hash<std::thread::id>{}(std::this_thread::get_id()) + timestamp);

    for (unsigned i = 0; i < pop_size; i++) {
        for (unsigned j = 0; j < dim; j++) {
            pop[i][j] = dist(gen);
        }
    }
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

void Ga::shade() {
    Ga::prev_min_fitness = Ga::min_fitness;

    for (int i = 0; i < pop_size; i++) {
        generate_trial_vector(trial_pop[i], i, trial_cr[i], trial_f[i]);
        trial_fitness[i] = function1(trial_pop[i]);
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

    rate_of_improvement = (Ga::prev_min_fitness - Ga::min_fitness) / Ga::prev_min_fitness;
    printf("Rate of improvement: %f\n", rate_of_improvement);
    if (rate_of_improvement < 0.001) {
        Ga::convergence_counter++;
    } else {
        Ga::convergence_counter = 0;
    }
}

void Ga::evolve(int generations, bool break_on_convergence) {
    int i;
    Ga::convergence_counter = 0;
    float convergence_threshold = 0.1;
    for (i=0; i<generations; i++) {
        compute_fitness();
        Ga::min_fitness_vector.push_back(min_fitness);

        printf("%i: %e\n", i, min_fitness);

        selection_roulette();

        crossover_uniform();

        mutation_normal_dist();

        if(break_on_convergence){
            if (convergence < convergence_threshold) {
                Ga::convergence_counter++;
                if (Ga::convergence_counter >= 10){
                    break;
                }
            } else {
                Ga::convergence_counter = 0;
            }
        }
    }

    printf("Generations: %d \n", i);
    printf("Best fitness: %e \n", Ga::best_fitness);
}

void Ga::evolve_shade(int generations) {
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


