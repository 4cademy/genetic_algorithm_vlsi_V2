//
// Created by Marcel Beyer on 22.08.2023.
//

#include "initialise.h"
#include <iostream>
#include <fstream>
#include <sstream>

auto* opt1 = new double[1000];
int i = 0;

void load_data_f1() {
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
