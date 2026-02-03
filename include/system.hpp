#pragma once
#include <omp.h>
#include <memory>
#include <vector>
#include <random>

class Solver;
class Engine;

enum SolverType {
    SOLVER_EULER,
    SOLVER_CAUCHY_EULER,
    SOLVER_POSITION_VERLET,
    SOLVER_VELOCITY_VERLET,
    SOLVER_RK4,
    SOLVER_RK45,
    SOLVER_PEFRL
};

enum EngineType {
    ENGINE_DIRECT_FORCE,
    ENGINE_BARNES_HUT,
    ENGINE_FAST_MULTIPOLE
};

class System {
    friend class Engine;
    friend class Solver;

    protected:
    std::vector<double> x, y, z;
    std::vector<double> vx, vy, vz;
    std::vector<double> ax, ay, az;
    std::vector<double> m;
    int n;

    std::unique_ptr<Solver> solver;
    std::unique_ptr<Engine> engine;

    public:
    virtual ~System();
    System(SolverType s, EngineType e, int n);
    double getBodyMass(int id);
    int getBodyCount();
    bool switchSolver(SolverType s);
    bool switchEngine(EngineType e);
    void update();
    void step(double h);
    void run(double t, double h);
};

class Cluster : public System {
    private:
    static thread_local std::mt19937 gen;
    static thread_local std::uniform_real_distribution<double> dist;

    double getX();

    public:
    Cluster(SolverType s, EngineType e, int n);
    void runStats(double t, double h);
};