#include "system.hpp"
#include "engine.hpp"
#include "solver.hpp"
#include <omp.h>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <memory>

System::~System() = default;

System::System(SolverType s, EngineType e, int n) : n(n) {
    x.resize(n);
    y.resize(n);
    z.resize(n);
    vx.resize(n);
    vy.resize(n);
    vz.resize(n);
    ax.resize(n);
    ay.resize(n);
    az.resize(n);
    m.resize(n);

    bool success;
    success = switchSolver(s);
    if (!success) throw std::runtime_error("Could not initialize solver for the system");
    success = switchEngine(e);
    if (!success) throw std::runtime_error("Could not initialize engine for the system");
}

double System::getBodyMass(int id) {
    return m[id];
}

int System::getBodyCount() {
    return n;
}

bool System::switchSolver(SolverType s) {
    switch(s) {
    case SOLVER_EULER:
        solver = std::make_unique<Euler>(this);
        break;
    case SOLVER_CAUCHY_EULER:
        solver = std::make_unique<CauchyEuler>(this);
        break;
    case SOLVER_POSITION_VERLET:
        solver = std::make_unique<PositionVerlet>(this);
        break;
    case SOLVER_VELOCITY_VERLET:
        solver = std::make_unique<VelocityVerlet>(this);
        break;
    case SOLVER_RK4:
        solver = std::make_unique<RK4>(this);
        break;
    case SOLVER_RK45:
        solver = std::make_unique<RK45>(this);
        break;
    case SOLVER_PEFRL:
        solver = std::make_unique<PEFRL>(this);
        break;
    default:
        return false;
    }

    return true;
}

bool System::switchEngine(EngineType e) {
    switch(e) {
    case ENGINE_DIRECT_FORCE:
        engine = std::make_unique<DirectForce>(this);
        break;
    case ENGINE_BARNES_HUT:
        engine = std::make_unique<BarnesHut>(this);
        break;
    case ENGINE_FAST_MULTIPOLE:
        engine = std::make_unique<FastMultipole>(this);
        break;
    default:
        return false;
    }

    return true;
}

void System::update() {
    engine->update();
}

void System::step(double h) {
    solver->step(h);
}

void System::run(double t, double h) {
    int steps = static_cast<int>(t / h);
    for (int i = 0; i < steps; i++) step(h);
}

double Cluster::getX() {
    return dist(gen);
}

thread_local std::mt19937 Cluster::gen;
thread_local std::uniform_real_distribution<double> Cluster::dist{0.0, 1.0};

Cluster::Cluster(SolverType s, EngineType e, int n) : System(s, e, n) {
    double mass = 1.0 / n;

    double svx = 0.0, svy = 0.0, svz = 0.0;

    #pragma omp parallel for schedule(dynamic, 64) reduction(+:svx, svy, svz)
    for (int i = 0; i < n; i++) {
        m[i] = mass;

        double x1 = getX();
        double r = 1.0 * std::sqrt(1.0 / (std::cbrt(x1 * x1) - 1.0));

        double x2 = getX(), x3 = getX();
        z[i] = (1.0 - 2.0 * x2) * r;
        double a1 = std::sqrt(r * r - z[i] * z[i]);
        double theta = 2.0 * M_PI * x3;
        y[i] = a1 * std::sin(theta);
        x[i] = a1 * std::cos(theta);

        double ve = std::sqrt(2.0) * std::pow(1.0 + r * r, -0.25);
        double g, q;
        do {
            q = getX();
            g = q * q * std::pow(1.0 - q * q, 3.5);
        } while(0.1 * getX() >= g);
        
        double v = ve * q;
        double x6 = getX(), x7 = getX();
        vz[i] = (1.0 - 2.0 * x6) * v;
        double a2 = std::sqrt(v * v - vz[i] * vz[i]);
        double phi = 2.0 * M_PI * x7;
        vy[i] = a2 * std::sin(phi);
        vx[i] = a2 * std::cos(phi);

        svx += vx[i];
        svy += vy[i];
        svz += vz[i];
    }

    svx /= static_cast<double>(n);
    svy /= static_cast<double>(n);
    svz /= static_cast<double>(n);

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        vx[i] -= svx;
        vy[i] -= svy;
        vz[i] -= svz;
    }
}

void Cluster::runStats(double t, double h) {
    int steps = static_cast<int>(t / h);
    int width = snprintf(nullptr, 0, "%d", steps);
    for (int i = 1; i <= steps; i++) {
        step(h);
        if (i % 100 == 0 || i == steps) {
            printf("\rStep %*d / %*d", width, i, width, steps);
            fflush(stdout);
        }
    }
    printf("\nFinished!\n");
}