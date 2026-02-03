#include "solver.hpp"
#include "system.hpp"
#include <omp.h>
#include <cmath>

Solver::Solver(System* s) : s(s) {
    x = s->x.data();
    y = s->y.data();
    z = s->z.data();
    vx = s->vx.data();
    vy = s->vy.data();
    vz = s->vz.data();
    ax = s->ax.data();
    ay = s->ay.data();
    az = s->az.data();
    n = s->n;
}

Euler::Euler(System* s) : Solver(s) {}

void Euler::step(double h) {
    s->update();

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        x[i] = std::fma(vx[i], h, x[i]);
        y[i] = std::fma(vy[i], h, y[i]);
        z[i] = std::fma(vz[i], h, z[i]);

        vx[i] = std::fma(ax[i], h, vx[i]);
        vy[i] = std::fma(ay[i], h, vy[i]);
        vz[i] = std::fma(az[i], h, vz[i]);
    }
}

CauchyEuler::CauchyEuler(System* s) : Solver(s) {}

void CauchyEuler::step(double h) {
    s->update();

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        vx[i] = std::fma(ax[i], h, vx[i]);
        vy[i] = std::fma(ay[i], h, vy[i]);
        vz[i] = std::fma(az[i], h, vz[i]);

        x[i] = std::fma(vx[i], h, x[i]);
        y[i] = std::fma(vy[i], h, y[i]);
        z[i] = std::fma(vz[i], h, z[i]);
    }
}

PositionVerlet::PositionVerlet(System* s) : Solver(s) {}

void PositionVerlet::step(double h) {
    double half = 0.5 * h;

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        x[i] = std::fma(vx[i], half, x[i]);
        y[i] = std::fma(vy[i], half, y[i]);
        z[i] = std::fma(vz[i], half, z[i]);
    }

    s->update();

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        vx[i] = std::fma(ax[i], h, vx[i]);
        vy[i] = std::fma(ay[i], h, vy[i]);
        vz[i] = std::fma(az[i], h, vz[i]);

        
        x[i] = std::fma(vx[i], half, x[i]);
        y[i] = std::fma(vy[i], half, y[i]);
        z[i] = std::fma(vz[i], half, z[i]);
    }
}

VelocityVerlet::VelocityVerlet(System* s) : Solver(s) {}

void VelocityVerlet::step(double) {
    
}

RK4::RK4(System* s) : Solver(s) {}

void RK4::step(double) {

}

RK45::RK45(System* s) : RK4(s) {}

void RK45::step(double) {

}

PEFRL::PEFRL(System* s) : Solver(s) {}

void PEFRL::step(double) {

}
