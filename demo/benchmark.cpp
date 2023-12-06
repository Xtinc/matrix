#include "benchmark/benchmark.h"
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

BENCHMARK(gemm_for_3M);
BENCHMARK(gemm_for_4M);
BENCHMARK_MAIN();
