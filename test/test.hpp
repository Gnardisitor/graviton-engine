#include "system.hpp"

class TestSystem : public System {
    private:
    static constexpr double masses[] = {1.989e30, 5.972e24};
    static constexpr double positions[] = {-0.000673, 0.000222, -0.000027, -0.178961, 0.967260, -0.000006};
    static constexpr double velocities[] = {0.000006, -0.000000, -0.000000, -0.172668, 0.968633, -0.000040};
    static constexpr double expectedPositions[] = {0.001151, 0.000822, -0.000059, -0.017198, -0.003189, -0.000001};
    static constexpr double expectedVelocities[] = {0.000004, 0.000003, -0.000000, -0.017202, -0.003101, 0.000001};

    static constexpr double g = 6.674e-11;
    static constexpr double au = 1.496e11;
    static constexpr double year = 365.0 * 86400.0;

    static constexpr double massSun = 1.989e30;
    static constexpr double massScale = (g * year * year) / (au * au * au * massSun * massSun);

    public:
    TestSystem(SolverType s, EngineType e);
    void resetVectors(int index);
    double getEnergy();
    bool runTest(double step, double delta);
};