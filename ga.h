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

    void clean() const;

    double** pop{};
    double* fitness{};
    double min_fitness{};
    double max_fitness{};
    double* opt1{};

private:
    unsigned dim{};
    unsigned pop_size{};
    double min_gene{};
    double max_gene{};

    double function1(const double* individual) const;
    void compute_fitness();
};


#endif //GENETIC_ALGORITHM_VLSI_V2_GA_H
