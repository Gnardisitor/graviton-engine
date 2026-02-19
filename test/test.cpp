#include "test.hpp"
#include "system.hpp"
#include <stdexcept>
#include <cmath>

TestSystem::TestSystem(SolverType s, EngineType e) : System(s, e, 2) {
    for (int i = 0; i < n; i++) {
        m[i] = masses[i] * massScale;
        x[i] = positions[i * 3 + 0];
        y[i] = positions[i * 3 + 1];
        z[i] = positions[i * 3 + 2];
        vx[i] = velocities[i * 3 + 0];
        vy[i] = velocities[i * 3 + 1];
        vz[i] = velocities[i * 3 + 2];
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

    if (std::fabs(finalEnergy - initialEnergy) > delta) return false;
    return true;
}