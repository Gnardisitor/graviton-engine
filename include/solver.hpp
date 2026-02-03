#pragma once

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
    // Add required variables and functions

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
    // Add required variables and functions
    
    public:
    PEFRL(System* s);
    void step(double h);
};