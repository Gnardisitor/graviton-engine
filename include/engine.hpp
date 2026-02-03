#pragma once
#include "system.hpp"

class Engine {
    protected:
    System* s;
    double *__restrict__ x, *__restrict__ y, *__restrict__ z;
    double *__restrict__ vx, *__restrict__ vy, *__restrict__ vz;
    double *__restrict__ ax, *__restrict__ ay, *__restrict__ az;
    double *__restrict__ m;
    int n;

    public:
    virtual ~Engine() = default;
    Engine(System* s);
    virtual void update() = 0;
};

class DirectForce : public Engine {
    public:
    DirectForce(System* s);
    void update();
};

class BarnesHut : public Engine {
    private:
    // Add required variables and functions

    public:
    BarnesHut(System * s);
    void update();
};

class FastMultipole : public Engine {
    private:
    // Add required variables and functions

    public:
    FastMultipole(System* s);
    void update();
};