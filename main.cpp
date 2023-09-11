//
// Created by Marcel Beyer on 22.08.2023.
//
#include <cstdio>
#include "ga.h"
#include <chrono>
#include <cstdlib>

using namespace std::chrono;

unsigned dim = 1000;
unsigned pop_size = 1'000;

int main() {
    auto start = high_resolution_clock::now();
    printf("%i\n", RAND_MAX);
    Ga* ga = new Ga(dim, pop_size, -100, 100);

    ga->evolve(1000);
    // free allocated memory
    ga->clean();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    printf("Time taken by function: %lld seconds", duration.count());
    return 0;
}