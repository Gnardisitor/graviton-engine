#include "main.hpp"
#include "system.hpp"
#include <omp.h>
#include <memory>
#include <chrono>
#include <iostream>

int main() {
    if (THREADS > 0) omp_set_num_threads(THREADS);

    std::cout << "Graviton test with " << COUNT << " bodies\n";

    auto start = std::chrono::high_resolution_clock::now();
    std::unique_ptr<Cluster> cluster = std::make_unique<Cluster>(SOLVER, ENGINE, COUNT);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Initialization took " << duration.count() << " s\n";

    start = std::chrono::high_resolution_clock::now();
    cluster->runStats(TIME, STEP);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "Simulation took " << duration.count() << " s\n";

    return 0;
}