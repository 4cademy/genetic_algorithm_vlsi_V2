//
// Created by Marcel Beyer on 22.08.2023.
//
#include <random>
#include "ga.h"
#include <chrono>
using namespace std::chrono;

unsigned dim = 1000;
unsigned pop_size = 1000;

int main() {
    auto start = high_resolution_clock::now();
    Ga* ga = new Ga(dim, pop_size, -100, 100);

    ga->evolve(1000);

    // free allocated memory
    ga->clean();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    printf("Time taken by function: %lld microseconds", duration.count());
    return 0;
}