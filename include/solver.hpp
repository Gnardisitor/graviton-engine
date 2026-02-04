#pragma once

#include <vector>
class System;

class Solver {
    protected:
    System* s;
    double *__restrict__ x, *__restrict__ y, *__restrict__ z;
    double *__restrict__ vx, *__restrict__ vy, *__restrict__ vz;
    double *__restrict__ ax, *__restrict__ ay, *__restrict__ az;
    int n;

    public:
    virtual ~Solver() = default;
    Solver(System* s);
    virtual void step(double h) = 0;
};

class Euler : public Solver {
    public:
    Euler(System* s);
    void step(double h);
};

class CauchyEuler : public Solver {
    public:
    CauchyEuler(System* s);
    void step(double h);
};

class PositionVerlet : public Solver {
    public:
    PositionVerlet(System* s);
    void step(double h);
};

class VelocityVerlet : public Solver {
    public:
    VelocityVerlet(System* s);
    void step(double h);
};

class RK4 : public Solver {
    protected:
    std::vector<double> ytx, yty, ytz, ytvx, ytvy, ytvz;
    std::vector<double> k1x, k1y, k1z, k1vx, k1vy, k1vz;
    std::vector<double> k2x, k2y, k2z, k2vx, k2vy, k2vz;
    std::vector<double> k3x, k3y, k3z, k3vx, k3vy, k3vz;
    std::vector<double> k4x, k4y, k4z, k4vx, k4vy, k4vz;

    void f(std::vector<double>& kx, std::vector<double>& ky, std::vector<double>& kz,
           std::vector<double>& kvx, std::vector<double>& kvy, std::vector<double>& kvz);

    public:
    RK4(System* s);
    void step(double h);
};

class RK45 : public RK4 {
    public:
    RK45(System* s);
    void step(double h);
};

class PEFRL : public Solver {
    private:
    static constexpr double XI = 0.1786178958448091;
    static constexpr double LAMBDA = -0.2123418310626054;
    static constexpr double CHI = -0.6626458266981849e-1;
    static constexpr double V1 = 1.0 - 2.0 * LAMBDA;
    static constexpr double V2 = 1.0 - 2.0 * (CHI + XI);
    
    public:
    PEFRL(System* s);
    void step(double h);
};