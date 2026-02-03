#include "engine.hpp"
#include <omp.h>
#include <cmath>

Engine::Engine(System *s) : s(s) {
    x = s->x.data();
    y = s->y.data();
    z = s->z.data();
    vx = s->vx.data();
    vy = s->vy.data();
    vz = s->vz.data();
    ax = s->ax.data();
    ay = s->ay.data();
    az = s->az.data();
    m = s->m.data();
    n = s->n;
}

DirectForce::DirectForce(System* s) : Engine(s) {}

void DirectForce::update() {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        double lax = 0.0, lay = 0.0, laz = 0.0;

        #pragma omp simd reduction(+:lax, lay, laz)
        for (int j = 0; j < i; j++) {
            double dx = x[j] - x[i];
            double dy = y[j] - y[i];
            double dz = z[j] - z[i];

            double r = 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz);
            double d = r * r * r;
            
            double f = m[j] * d;

            lax = std::fma(dx, f, lax);
            lay = std::fma(dy, f, lay);
            laz = std::fma(dz, f, laz);
        }

        #pragma omp simd reduction(+:lax, lay, laz)
        for (int j = i + 1; j < n; j++) {
            double dx = x[j] - x[i];
            double dy = y[j] - y[i];
            double dz = z[j] - z[i];

            double r = 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz);
            double d = r * r * r;
            
            double f = m[j] * d;

            lax = std::fma(dx, f, lax);
            lay = std::fma(dy, f, lay);
            laz = std::fma(dz, f, laz);
        }

        ax[i] = lax;
        ay[i] = lay;
        az[i] = laz;
    }
}

BarnesHut::BarnesHut(System* s) : Engine(s) {}

void BarnesHut::update() {

}

FastMultipole::FastMultipole(System* s) : Engine(s) {}

void FastMultipole::update() {
    
}