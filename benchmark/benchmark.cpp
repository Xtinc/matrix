#include "benchmark/benchmark.h"
#include "../demo/robotics.hpp"
#include "liegroup.hpp"

using namespace ppx;

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
}

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
    // std::cout << "singualr q:" << q << "\n";
}

BENCHMARK(gemm_for_3M);
BENCHMARK(gemm_for_4M);
BENCHMARK(codo_perf);
BENCHMARK_MAIN();
