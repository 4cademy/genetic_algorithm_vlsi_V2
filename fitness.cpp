//
// Created by Marcel Beyer on 04.09.2023.
//
#include "fitness.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <omp.h>
#include <cmath>

extern unsigned dim;
double* opt1;

void load_data_f1() {
    int i = 0;
    std::string path = "../cdatafiles/F1-xopt.txt";
    std::ifstream file(path);

    if (file.is_open()) {
        std::string line;
        while (getline(file, line)) {
            opt1[i++] = stod(line);
        }
        file.close();
    }
    else {
        std::cout << "Cannot open the datafile '" << path << "'" << std::endl;
    }
}

Fitness::Fitness(){
    opt1 = new double[dim];
    load_data_f1();
}

void Fitness::debug_print() {
    printf("dim: %d\n", dim);
    for(int i = 0; i < dim; i++) {
        printf("opt1[%d]: %f\n", i, opt1[i]);
    }
}

double Fitness::function1(const double* x) {
    double result = 0;
    auto* z = new double[dim];

    int sign;
    double hat;
    double c1;
    double c2;

// # pragma omp parallel for default(none) shared(x, z, dim, opt1) private(sign, hat, c1, c2) reduction(+:result)
    for(unsigned i = 0; i < dim; i++) {
        z[i] = x[i] - opt1[i];
        // Transformation
        if (z[i] > 0) {
            sign = 1;
            hat = log(fabs(z[i]));
            c1 = 10;
            c2 = 7.9;
            z[i] = sign*exp(hat+0.049*(sin(c1*hat)+sin(c2*hat)));
        } else if (z[i] == 0) {
            z[i] = 0;
        } else {
            sign = -1;
            hat = log(fabs(z[i]));
            c1 = 5.5;
            c2 = 3.1;
            z[i] = sign*exp(hat+0.049*(sin(c1*hat)+sin(c2*hat)));
        }
        result += pow(1.0e6,  i/((double)(dim - 1)) ) * z[i] * z[i];
    }

    delete[] z;
    return result;
}