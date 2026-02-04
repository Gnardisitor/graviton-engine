#include "solver.hpp"
#include "system.hpp"
#include <omp.h>
#include <cmath>
#include <stdexcept>

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

void VelocityVerlet::step(double h) {
    double half = 0.5 * h;

    s->update();

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        vx[i] = std::fma(ax[i], half, vx[i]);
        vy[i] = std::fma(ay[i], half, vy[i]);
        vz[i] = std::fma(az[i], half, vz[i]);
    }

    s->update();

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        x[i] = std::fma(vx[i], h, x[i]);
        y[i] = std::fma(vy[i], h, y[i]);
        z[i] = std::fma(vz[i], h, z[i]);

        vx[i] = std::fma(ax[i], half, vx[i]);
        vy[i] = std::fma(ay[i], half, vy[i]);
        vz[i] = std::fma(az[i], half, vz[i]);
    }
}

RK4::RK4(System* s) : Solver(s) {
    ytx.resize(n);
    yty.resize(n);
    ytz.resize(n);
    ytvx.resize(n);
    ytvy.resize(n);
    ytvz.resize(n);
    k1x.resize(n);
    k1y.resize(n);
    k1z.resize(n);
    k1vx.resize(n);
    k1vy.resize(n);
    k1vz.resize(n);
    k2x.resize(n);
    k2y.resize(n);
    k2z.resize(n);
    k2vx.resize(n);
    k2vy.resize(n);
    k2vz.resize(n);
    k3x.resize(n);
    k3y.resize(n);
    k3z.resize(n);
    k3vx.resize(n);
    k3vy.resize(n);
    k3vz.resize(n);
    k4x.resize(n);
    k4y.resize(n);
    k4z.resize(n);
    k4vx.resize(n);
    k4vy.resize(n);
    k4vz.resize(n);
}

void RK4::f(std::vector<double>& kx, std::vector<double>& ky, std::vector<double>& kz,
            std::vector<double>& kvx, std::vector<double>& kvy, std::vector<double>& kvz) {
    s->update();
    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        kx[i] = vx[i];
        ky[i] = vy[i];
        kz[i] = vz[i];
        kvx[i] = ax[i];
        kvy[i] = ay[i];
        kvz[i] = az[i];
    }
}

void RK4::step(double h) {
    double half = 0.5 * h;
    double sixth = h / 6.0;
    
    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        ytx[i] = x[i];
        yty[i] = y[i];
        ytz[i] = z[i];
        ytvx[i] = vx[i];
        ytvy[i] = vy[i];
        ytvz[i] = vz[i];
    }
    
    f(k1x, k1y, k1z, k1vx, k1vy, k1vz);

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        x[i] = std::fma(k1x[i], half, ytx[i]);
        y[i] = std::fma(k1y[i], half, yty[i]);
        z[i] = std::fma(k1z[i], half, ytz[i]);
        vx[i] = std::fma(k1vx[i], half, ytvx[i]);
        vy[i] = std::fma(k1vy[i], half, ytvy[i]);
        vz[i] = std::fma(k1vz[i], half, ytvz[i]);
    }
    
    f(k2x, k2y, k2z, k2vx, k2vy, k2vz);

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        x[i] = std::fma(k2x[i], half, ytx[i]);
        y[i] = std::fma(k2y[i], half, yty[i]);
        z[i] = std::fma(k2z[i], half, ytz[i]);
        vx[i] = std::fma(k2vx[i], half, ytvx[i]);
        vy[i] = std::fma(k2vy[i], half, ytvy[i]);
        vz[i] = std::fma(k2vz[i], half, ytvz[i]);
    }
    
    f(k3x, k3y, k3z, k3vx, k3vy, k3vz);

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        x[i] = std::fma(k3x[i], h, ytx[i]);
        y[i] = std::fma(k3y[i], h, yty[i]);
        z[i] = std::fma(k3z[i], h, ytz[i]);
        vx[i] = std::fma(k3vx[i], h, ytvx[i]);
        vy[i] = std::fma(k3vy[i], h, ytvy[i]);
        vz[i] = std::fma(k3vz[i], h, ytvz[i]);
    }
    
    f(k4x, k4y, k4z, k4vx, k4vy, k4vz);

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        double sum_kx = k1x[i] + 2.0*k2x[i] + 2.0*k3x[i] + k4x[i];
        double sum_ky = k1y[i] + 2.0*k2y[i] + 2.0*k3y[i] + k4y[i];
        double sum_kz = k1z[i] + 2.0*k2z[i] + 2.0*k3z[i] + k4z[i];
        double sum_kvx = k1vx[i] + 2.0*k2vx[i] + 2.0*k3vx[i] + k4vx[i];
        double sum_kvy = k1vy[i] + 2.0*k2vy[i] + 2.0*k3vy[i] + k4vy[i];
        double sum_kvz = k1vz[i] + 2.0*k2vz[i] + 2.0*k3vz[i] + k4vz[i];
        
        x[i] = std::fma(sixth, sum_kx, ytx[i]);
        y[i] = std::fma(sixth, sum_ky, yty[i]);
        z[i] = std::fma(sixth, sum_kz, ytz[i]);
        vx[i] = std::fma(sixth, sum_kvx, ytvx[i]);
        vy[i] = std::fma(sixth, sum_kvy, ytvy[i]);
        vz[i] = std::fma(sixth, sum_kvz, ytvz[i]);
    }
}

RK45::RK45(System* s) : RK4(s) {}

void RK45::step(double) {
    throw std::runtime_error("RK45 is not implemented");
}

PEFRL::PEFRL(System* s) : Solver(s) {}

void PEFRL::step(double h) {
    double half = 0.5 * h;
    double v1 = XI * h;
    double v2 = V1 * half;
    double v3 = LAMBDA * h;
    double v4 = V2 * half;
    
    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        x[i] = std::fma(vx[i], v1, x[i]);
        y[i] = std::fma(vy[i], v1, y[i]);
        z[i] = std::fma(vz[i], v1, z[i]);
    }

    s->update();

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        vx[i] = std::fma(ax[i], v2, vx[i]);
        vy[i] = std::fma(ay[i], v2, vy[i]);
        vz[i] = std::fma(az[i], v2, vz[i]);

        x[i] = std::fma(vx[i], v1, x[i]);
        y[i] = std::fma(vy[i], v1, y[i]);
        z[i] = std::fma(vz[i], v1, z[i]);
    }

    s->update();

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        vx[i] = std::fma(ax[i], v3, vx[i]);
        vy[i] = std::fma(ay[i], v3, vy[i]);
        vz[i] = std::fma(az[i], v3, vz[i]);

        x[i] = std::fma(vx[i], v4, x[i]);
        y[i] = std::fma(vy[i], v4, y[i]);
        z[i] = std::fma(vz[i], v4, z[i]);
    }

    s->update();

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        vx[i] = std::fma(ax[i], v3, vx[i]);
        vy[i] = std::fma(ay[i], v3, vy[i]);
        vz[i] = std::fma(az[i], v3, vz[i]);

        x[i] = std::fma(vx[i], v1, x[i]);
        y[i] = std::fma(vy[i], v1, y[i]);
        z[i] = std::fma(vz[i], v1, z[i]);
    }

    s->update();

    #pragma omp parallel for simd schedule(static)
    for (int i = 0; i < n; i++) {
        vx[i] = std::fma(ax[i], v2, vx[i]);
        vy[i] = std::fma(ay[i], v2, vy[i]);
        vz[i] = std::fma(az[i], v2, vz[i]);

        x[i] = std::fma(vx[i], v1, x[i]);
        y[i] = std::fma(vy[i], v1, y[i]);
        z[i] = std::fma(vz[i], v1, z[i]);
    }
}
