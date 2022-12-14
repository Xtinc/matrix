#include "Eigen/Dense"
#include "modern_robotics.h"
#include <benchmark/benchmark.h>
#include "LieGroup.hpp"
#include <random>

using namespace ppx;

static void BM_MatrixMul(benchmark::State &state)
{
    std::default_random_engine e;
    std::uniform_real_distribution<> u(-1.0e6, 1.0e6);
    Matrix<20, 12> M{};
    for (size_t i = 0; i < 20; i++)
    {
        for (size_t j = 0; i < 12; i++)
        {
            double t = u(e);
            M(i, j) = t;
        }
    }
    auto b = M.T();
    Matrix<20, 20> M2{};
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
    Eigen::MatrixXd M2(20, 12);
    for (size_t i = 0; i < 20; i++)
    {
        for (size_t j = 0; i < 12; i++)
        {
            double t = u(e);
            M2(i, j) = t;
        }
    }

    auto b = M2.transpose();
    Eigen::MatrixXd M3(20, 20);
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

static void BM_MatrixExpr(benchmark::State &state)
{
    Matrix<10, 10> x =
        std::array<double, 100>{
            -729046.0, 451678.0, 595857.0, 594560.0, -748207.0, 411549.0, 642492.0, 977043.0, -864809.0, -191583.0, 670017., 962219., -277412, -366899, -579582, -994363, 641682, -335103, 587195, -294475, 937736, -780276, -576151, 744858, -897567, 421408, 880148, -400337, 189007, 185648 - 557932, 596212, 362719, -701772, -927117, 287922, -174667, -972922, 465597, -287310 - 383666, -405941, -202523, 988137, -182538, -87934.4, -153670, -565524, 390466, 929933, 94441.2, -990433, 481294, 643807, -84021.7, 547834, 161913, 814729, 359640, -691123 - 623236, -775071, -50482.6, -749634, -24862.1, 147509, -683885, 696936, -215359, -210184, 985763, 279527, -155825, 527500, 587950, 753515, 523462, 910035, 123115, -225408, 992923, 756861, -652270, -18821.9, 841750, 616351, -539688, 557795, -583864, 453909, 935390, 7325.36, -396174, 327211, 615062, -964452, 619469, 974919, 54742.9, -222860};
    auto y = x.T();
    Matrix<10, 10> z{};
    for (auto _ : state)
    {
        z = x * 2 + y * 2 + 3.3 + x * y + x + y + x + y;
    }
}

static void BM_MatrixEigenExpr(benchmark::State &state)
{
    Eigen::MatrixXd x(10, 10);
    x << -729046, 451678, 595857, 594560, -748207, 411549, 642492, 977043, -864809, -191583, 670017, 962219, -277412, -366899, -579582, -994363, 641682, -335103, 587195, -294475, 937736, -780276, -576151, 744858, -897567, 421408, 880148, -400337, 189007, 185648 - 557932, 596212, 362719, -701772, -927117, 287922, -174667, -972922, 465597, -287310 - 383666, -405941, -202523, 988137, -182538, -87934.4, -153670, -565524, 390466, 929933, 94441.2, -990433, 481294, 643807, -84021.7, 547834, 161913, 814729, 359640, -691123 - 623236, -775071, -50482.6, -749634, -24862.1, 147509, -683885, 696936, -215359, -210184, 985763, 279527, -155825, 527500, 587950, 753515, 523462, 910035, 123115, -225408, 992923, 756861, -652270, -18821.9, 841750, 616351, -539688, 557795, -583864, 453909, 935390, 7325.36, -396174, 327211, 615062, -964452, 619469, 974919, 54742.9, -222860;
    Eigen::MatrixXd y = x.transpose();
    Eigen::MatrixXd z(10, 10);
    Eigen::MatrixXd s3(10, 10);
    s3.fill(3.3);
    for (auto _ : state)
    {
        z = x * 2 + y * 2 + s3 + x * y + x + y + x + y;
    }
}

static void BM_MatrixLog(benchmark::State &state)
{
    SE3 SE3mat = {1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 3, 1};
    se3 tmp{};
    for (auto _ : state)
    {
        tmp = SE3mat.log();
    }
}

static void BM_MatrixExp(benchmark::State &state)
{
    Matrix<4, 4> se3mat = {0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 1.5708, 0.0,
                           0.0, -1.5708, 0.0, 0.0,
                           0.0, 2.3562, 2.3562, 0.0};
    SE3 tmp{};
    for (auto _ : state)
    {
        tmp = vee(se3mat).exp();
    }
}

static void BM_MatrixEigenLog(benchmark::State &state)
{
    Eigen::MatrixXd result(4, 4);
    Eigen::MatrixXd tmp{};
    result << 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 3, 0, 0, 0, 1;
    for (auto _ : state)
    {
        tmp = mr::MatrixLog6(result);
    }
}

static void BM_MatrixEigenExp(benchmark::State &state)
{
    Eigen::MatrixXd result(4, 4);
    Eigen::MatrixXd tmp{};
    result << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.5708, 2.3562, 0.0, 1.5708, 0.0, 2.3562, 0.0, 0.0, 0.0, 0.0;
    for (auto _ : state)
    {
        tmp = mr::MatrixExp6(result);
    }
}

static void BM_MatrixQR(benchmark::State &state)
{
    Matrix<5, 6> u{1233.0, 415.0, 87.7, 11.6, 243.0,
                   997.0, -122.0, 35.4, 889.0, 111.1,
                   -442.0, -0.987, 355.0, -346.0, 3419.0,
                   235.0, 98.87, -827.0, 876.0, 34.0,
                   222.0, -87.8, 546.0, -101.0, 122.1,
                   -86.0, 999.0, 65.2, 902.0, 54.2};
    Matrix<6, 5> result;
    Matrix<5, 1> c;
    Matrix<5, 1> d;
    bool sing;
    auto uT = u.T();
    for (auto _ : state)
    {
        result = qrdcmp(uT, c, d, sing);
    }
}

static void BM_MatrixEigenQR(benchmark::State &state)
{
    Eigen::MatrixXd u(6, 5);
    u << 1233.0, 415.0, 87.7, 11.6, 243.0,
        997.0, -122.0, 35.4, 889.0, 111.1,
        -442.0, -0.987, 355.0, -346.0, 3419.0,
        235.0, 98.87, -827.0, 876.0, 34.0,
        222.0, -87.8, 546.0, -101.0, 122.1,
        -86.0, 999.0, 65.2, 902.0, 54.2;
    auto dd = Eigen::HouseholderQR<Eigen::MatrixXd>(u);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(u.rows(), u.cols());
    for (auto _ : state)
    {
        dd.compute(u);
    }
}

static void BM_MatrixSVD(benchmark::State &state)
{
    Matrix<5, 6> u{1233.0, 415.0, 87.7, 11.6, 243.0,
                   997.0, -122.0, 35.4, 889.0, 111.1,
                   -442.0, -0.987, 355.0, -346.0, 3419.0,
                   235.0, 98.87, -827.0, 876.0, 34.0,
                   222.0, -87.8, 546.0, -101.0, 122.1,
                   -86.0, 999.0, 65.2, 902.0, 54.2};
    Matrix<5, 6> result;
    Matrix<6, 1> e;
    Matrix<6, 6> v;
    bool sing = false;
    for (auto _ : state)
    {
        result = svdcmp(u, e, v, sing);
        auto sss = result * v;
    }
}

static void BM_MatrixEigenSVD(benchmark::State &state)
{
    Eigen::MatrixXd u(5, 6);
    u << 1233, 997, -442, 235, 222, -86,
        415, -122, -0.987, 98.87, -87.8, 999,
        87.7, 35.4, 355, -827, 546, 65.2,
        11.6, 889, -346, 876, -101, 902,
        243, 111.1, 3419, 34, 122.1, 54.2;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(u, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd U(5, 6), S, V(6, 6);
    for (auto _ : state)
    {
        U = svd.matrixU();
        S = svd.singularValues();
        V = svd.matrixV();
        Eigen::MatrixXd RES = U * S * V;
    }
}

// Register the function as a benchmark
BENCHMARK(BM_MatrixMul);
BENCHMARK(BM_MatrixMulEigen);
BENCHMARK(BM_MatrixInv);
BENCHMARK(BM_MatrixEigenInv);
BENCHMARK(BM_MatrixExpr);
BENCHMARK(BM_MatrixEigenExpr);
BENCHMARK(BM_MatrixLog);
BENCHMARK(BM_MatrixEigenLog);
BENCHMARK(BM_MatrixExp);
BENCHMARK(BM_MatrixEigenExp);
BENCHMARK(BM_MatrixQR);
BENCHMARK(BM_MatrixEigenQR);
BENCHMARK(BM_MatrixSVD);
BENCHMARK(BM_MatrixEigenSVD);
// Run the benchmark
BENCHMARK_MAIN();