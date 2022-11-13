#include "Eigen/Dense"
#include <benchmark/benchmark.h>
#include "matrix.hpp"
#include <random>

static void BM_MatrixMul(benchmark::State &state)
{
    std::default_random_engine e;
    std::uniform_real_distribution<> u(-1.0e6, 1.0e6);
    Matrix<500, 15> M{};
    for (size_t i = 0; i < 500; i++)
    {
        for (size_t j = 0; i < 15; i++)
        {
            double t = u(e);
            M(i, j) = t;
        }
    }
    auto b = M.T();
    Matrix<500, 500> M2{};
    for (auto _ : state)
    {
        M2 = M * b;
    }
}

static void BM_MatrixInv(benchmark::State &state)
{
    Matrix<10, 10> M =
        std::array<double, 100>{
            -729046.0, 451678.0, 595857.0, 594560.0, -748207.0, 411549.0, 642492.0, 977043.0, -864809.0, -191583.0, 670017., 962219., -277412, -366899, -579582, -994363, 641682, -335103, 587195, -294475, 937736, -780276, -576151, 744858, -897567, 421408, 880148, -400337, 189007, 185648 - 557932, 596212, 362719, -701772, -927117, 287922, -174667, -972922, 465597, -287310 - 383666, -405941, -202523, 988137, -182538, -87934.4, -153670, -565524, 390466, 929933, 94441.2, -990433, 481294, 643807, -84021.7, 547834, 161913, 814729, 359640, -691123 - 623236, -775071, -50482.6, -749634, -24862.1, 147509, -683885, 696936, -215359, -210184, 985763, 279527, -155825, 527500, 587950, 753515, 523462, 910035, 123115, -225408, 992923, 756861, -652270, -18821.9, 841750, 616351, -539688, 557795, -583864, 453909, 935390, 7325.36, -396174, 327211, 615062, -964452, 619469, 974919, 54742.9, -222860};
    M = M.T();
    Matrix<10, 10> M2{};
    for (auto _ : state)
    {
        M2 = inverse(M);
    }
}

static void BM_MatrixMulEigen(benchmark::State &state)
{
    std::default_random_engine e;
    std::uniform_real_distribution<> u(-1.0e6, 1.0e6);
    Eigen::MatrixXd M2(500, 15);
    for (size_t i = 0; i < 500; i++)
    {
        for (size_t j = 0; i < 15; i++)
        {
            double t = u(e);
            M2(i, j) = t;
        }
    }

    auto b = M2.transpose();
    Eigen::MatrixXd M3(500, 500);
    for (auto _ : state)
    {
        M3 = M2 * b;
    }
}

static void BM_MatrixEigenInv(benchmark::State &state)
{
    Eigen::MatrixXd M2(10, 10);
    M2 << -729046, 451678, 595857, 594560, -748207, 411549, 642492, 977043, -864809, -191583, 670017, 962219, -277412, -366899, -579582, -994363, 641682, -335103, 587195, -294475, 937736, -780276, -576151, 744858, -897567, 421408, 880148, -400337, 189007, 185648 - 557932, 596212, 362719, -701772, -927117, 287922, -174667, -972922, 465597, -287310 - 383666, -405941, -202523, 988137, -182538, -87934.4, -153670, -565524, 390466, 929933, 94441.2, -990433, 481294, 643807, -84021.7, 547834, 161913, 814729, 359640, -691123 - 623236, -775071, -50482.6, -749634, -24862.1, 147509, -683885, 696936, -215359, -210184, 985763, 279527, -155825, 527500, 587950, 753515, 523462, 910035, 123115, -225408, 992923, 756861, -652270, -18821.9, 841750, 616351, -539688, 557795, -583864, 453909, 935390, 7325.36, -396174, 327211, 615062, -964452, 619469, 974919, 54742.9, -222860;
    Eigen::MatrixXd M;
    for (auto _ : state)
    {
        M = M2.inverse();
    }
}

// Register the function as a benchmark
BENCHMARK(BM_MatrixMul);
BENCHMARK(BM_MatrixMulEigen);
BENCHMARK(BM_MatrixInv);
BENCHMARK(BM_MatrixEigenInv);
// Run the benchmark
BENCHMARK_MAIN();