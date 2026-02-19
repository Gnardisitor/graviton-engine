#include "test.hpp"
#include "system.hpp"
#include <stdexcept>
#include <cmath>
#include <omp.h>

TestSystem::TestSystem(SolverType s, EngineType e) : System(s, e, testBodies) {
    for (int i = 0; i < n; i++) {
        m[i] = masses[i] * massScale;
        x[i] = positions[i * 3 + 0];
        y[i] = positions[i * 3 + 1];
        z[i] = positions[i * 3 + 2];
        vx[i] = velocities[i * 3 + 0] * year;  // Convert AU/day to AU/year
        vy[i] = velocities[i * 3 + 1] * year;
        vz[i] = velocities[i * 3 + 2] * year;
    }

    double tm = 0.0;
    double tvx = 0.0, tvy = 0.0, tvz = 0.0;

    for (int i = 0; i < n; i++) {
        tm += m[i];
        tvx += m[i] * vx[i];
        tvy += m[i] * vy[i];
        tvz += m[i] * vz[i];
    }

    tvx /= tm;
    tvy /= tm;
    tvz /= tm;

    for (int i = 0; i < n; i++) {
        vx[i] -= tvx;
        vy[i] -= tvy;
        vz[i] -= tvz;
    }
}

void TestSystem::resetVectors(int index) {
    if (index == 0 || index == 1) {
        x[index] = positions[index * 3 + 0];
        y[index] = positions[index * 3 + 1];
        z[index] = positions[index * 3 + 2];
        vx[index] = velocities[index * 3 + 0];
        vy[index] = velocities[index * 3 + 1];
        vz[index] = velocities[index * 3 + 2];
    }
    else throw std::runtime_error("Only two bodies exist in the test system");
}

double TestSystem::getEnergy() {
    double energy = 0.0;

    for (int i = 0; i < n; i++) {
        energy += 0.5 * m[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = x[j] - x[i];
            double dy = y[j] - y[i];
            double dz = z[j] - z[i];

            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            energy -= (m[i] * m[j]) / r;
        }
    }

    return energy;
}

bool TestSystem::runTest(double h, double delta) {
    int steps = static_cast<int>(1.0 / h);

    double initialEnergy = getEnergy();
    for (int i = 0; i < steps; i++) step(h);
    double finalEnergy = getEnergy();

    double energyChange = std::fabs(finalEnergy - initialEnergy) / std::fabs(initialEnergy);
    if (energyChange > delta) return false;

    for (int i = 0; i < n; i++) {
        double dx = std::fabs(x[i] - expectedPositions[i * 3 + 0]);
        double dy = std::fabs(y[i] - expectedPositions[i * 3 + 1]);
        double dz = std::fabs(z[i] - expectedPositions[i * 3 + 2]);
        if (dx > delta || dy > delta || dz > delta) return false;
    }

    for (int i = 0; i < n; i++) {
        double dvx = std::fabs(vx[i] - expectedVelocities[i * 3 + 0] * year);
        double dvy = std::fabs(vy[i] - expectedVelocities[i * 3 + 1] * year);
        double dvz = std::fabs(vz[i] - expectedVelocities[i * 3 + 2] * year);
        if (dvx > delta || dvy > delta || dvz > delta) return false;
    }

    return true;
}

bool testSolver(SolverType s) {
    TestSystem test = TestSystem(s, ENGINE_DIRECT_FORCE);
    return test.runTest(STEP, DELTA);
}

int main(int argc, char* argv[]) {
    omp_set_num_threads(1);

    if (argc != 2) return 1;
    int solverIndex = std::stoi(argv[1]);
    if (solverIndex < 0 || solverIndex > 6) return 1;
    SolverType solverType = static_cast<SolverType>(solverIndex);

    if (solverType == SOLVER_RK45) return 1;
    bool passed = testSolver(solverType);
    
    if (passed) return 0;
    return 1;
}