#include "main.hpp"
#include "system.hpp"
#include <memory>
#include <iostream>

int main() {
    std::cout << "Graviton test with " << COUNT << " bodies\n";
    std::unique_ptr<Cluster> cluster = std::make_unique<Cluster>(SOLVER, ENGINE, COUNT);
    cluster->runStats(TIME, STEP);
    return 0;
}