#include "benchmark/benchmark.h"
#include "../demo/robotics.hpp"
#include "liegroup.hpp"

using namespace ppx;

static void gemm_for_3M(benchmark::State &state)
{
    MatrixS<3, 3> A;
    MatrixS<3, 3> B;
    random(A);
    random(B);
    MatrixS<3, 3> C;
    for (auto _ : state)
    {
        C = B * A;
    }
    benchmark::DoNotOptimize(C);
}

static void gemm_for_4M(benchmark::State &state)
{
    MatrixS<4, 4> A;
    MatrixS<4, 4> B;
    random(A);
    random(B);
    MatrixS<4, 4> C;
    for (auto _ : state)
    {
        C = B * A;
    }
    benchmark::DoNotOptimize(C);
}

static void gemm_for_5M(benchmark::State &state)
{
    MatrixS<5, 4> A;
    MatrixS<4, 5> B;
    random(A);
    random(B);
    MatrixS<5, 5> C;
    for (auto _ : state)
    {
        C = A * B;
    }
    benchmark::DoNotOptimize(C);
}

static void gemm_for_6M(benchmark::State &state)
{
    MatrixS<6, 4> A;
    MatrixS<4, 6> B;
    random(A);
    random(B);
    MatrixS<4, 4> C;
    for (auto _ : state)
    {
        C = B * A;
    }
    benchmark::DoNotOptimize(C);
}

static void gemm_for_9M(benchmark::State &state)
{
    MatrixS<9, 4> A;
    MatrixS<4, 7> B;
    random(A);
    random(B);
    MatrixS<9, 7> C;
    for (auto _ : state)
    {
        C = A * B;
    }
    benchmark::DoNotOptimize(C);
}

static void sum_for_5M(benchmark::State &state)
{
    MatrixS<5, 10> A;
    random(A);
    double s{};
    for (auto _ : state)
    {
        s = sum(A.data(), 50);
    }
    benchmark::DoNotOptimize(s);
}

static void sum_for_9M(benchmark::State &state)
{
    MatrixS<25, 10> A;
    random(A);
    double s{};
    for (auto _ : state)
    {
        s = sum(A.data(), 250);
    }
    benchmark::DoNotOptimize(s);
}

static void norm2_for_5M(benchmark::State &state)
{
    MatrixS<50, 1> A;
    random(A);
    double s{};
    for (auto _ : state)
    {
        s = norm2(A);
    }
    benchmark::DoNotOptimize(s);
}

static void norm2_for_9M(benchmark::State &state)
{
    MatrixS<250, 1> A;
    random(A);
    double s{};
    for (auto _ : state)
    {
        s = norm2(A);
    }
    benchmark::DoNotOptimize(s);
}

static void linsolve_LU(benchmark::State &state)
{
    MatrixS<16, 1> b;
    MatrixS<16, 16> A;
    random(A);
    b.fill(1);
    EqnResult<16> x;
    for (auto _ : state)
    {
        x = linsolve<Factorization::LU>(A, A * b);
    }
    benchmark::DoNotOptimize(x);
}

static void linsolve_QR(benchmark::State &state)
{
    MatrixS<16, 1> b;
    MatrixS<16, 16> A;
    random(A);
    b.fill(1);
    EqnResult<16> x;
    for (auto _ : state)
    {
        x = linsolve<Factorization::QR>(A, A * b);
    }
    benchmark::DoNotOptimize(x);
}

static void linsolve_SVD(benchmark::State &state)
{
    MatrixS<16, 1> b;
    MatrixS<16, 16> A;
    random(A);
    b.fill(1);
    EqnResult<16> x;
    for (auto _ : state)
    {
        x = linsolve<Factorization::SVD>(A, A * b);
    }
    benchmark::DoNotOptimize(x);
}

static void codo_perf(benchmark::State &state)
{
    kinematics<6> UR5;
    SE3 F6{-1.0, 0.0, 0.0, 0.0,
           0.0, 0.0, 1.0, 0.0,
           0.0, 1.0, 0.0, 0.0,
           0.817, 0.191, -0.006, 1.0};
    UR5.Joint<0>() = {"R1", se3{0, 0, 1, 0, 0, 0, 0}, SE3(), -PI, PI};
    UR5.Joint<1>() = {"R2", se3{0.0, 1.0, 0.0, -0.089, 0.0, 0.0}, SE3{}, -PI, PI};
    UR5.Joint<2>() = {"R3", se3{0.0, 1.0, 0.0, -0.089, 0.0, 0.425}, SE3{}, -PI, PI};
    UR5.Joint<3>() = {"R4", se3{0.0, 1.0, 0.0, -0.089, 0.0, 0.817}, SE3{}, -PI, PI};
    UR5.Joint<4>() = {"R5", se3{0.0, 0.0, -1.0, -0.109, 0.817, 0.0}, SE3{}, -PI, PI};
    UR5.Joint<5>() = {"R6", se3{0.0, 0.0, 0.0, 0.006, 0.0, 0.817}, F6, -PI, PI};
    kinematics<6>::Q q{-1.89406, 0.205429, 2.17377, -2.37914, 0.022822, 1.95129};
    // kinematics<6>::Q q{2.36254, -1.97934, -0.231841, -0.581726, 1.22112, 2.1659};

    SE3 TargetPose = UR5.forwardSpace(q);
    for (auto _ : state)
    {
        q = UR5.inverseSpace(TargetPose, q - 0.05);
    }
    benchmark::DoNotOptimize(q);
}

BENCHMARK(gemm_for_3M);
BENCHMARK(gemm_for_4M);
BENCHMARK(gemm_for_5M);
BENCHMARK(gemm_for_6M);
BENCHMARK(gemm_for_9M);
BENCHMARK(sum_for_5M);
BENCHMARK(sum_for_9M);
BENCHMARK(norm2_for_5M);
BENCHMARK(norm2_for_9M);
BENCHMARK(linsolve_LU);
BENCHMARK(linsolve_QR);
BENCHMARK(linsolve_SVD);
BENCHMARK(codo_perf);
BENCHMARK_MAIN();
