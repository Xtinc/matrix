#include "gtest/gtest.h"
#include "liegroup.hpp"

using namespace ppx;

class NonLinOpt_TestCase : public ::testing::Test
{
public:
    NonLinOpt_TestCase() = default;
};

TEST_F(NonLinOpt_TestCase, OneDimension_OPT)
{
    auto f1 = [](double x)
    {
        double y = 0.0;
        for (int k = -10; k < 11; k++)
        {
            y += (k + 1) * (k + 1) * cos(k * x) * exp(-1 * k * k / 2.0);
        }
        return y;
    };

    auto x1 = fminbnd<Optimization::GoldenSearch>(f1, 1.0, 3.0);
    auto x2 = fminbnd<Optimization::QuadraticSearch>(f1, 1.0, 3.0);
    auto x3 = fminbnd<Optimization::BrentSearch>(f1, 1.0, 3.0);
    EXPECT_NEAR(x1.x, 2.006062180123711, ppx::EPS_SP * 1e2);
    EXPECT_NEAR(x2.x, 2.006062180123711, ppx::EPS_SP * 1e2);
    EXPECT_NEAR(x3.x, 2.006062180123711, ppx::EPS_SP * 1e2);

    auto f2 = [](double x)
    {
        return sin(x);
    };

    x1 = fminbnd<Optimization::GoldenSearch>(f2, 0.0, 2.0 * PI);
    x2 = fminbnd<Optimization::QuadraticSearch>(f2, 0.0, 2.0 * PI);
    x3 = fminbnd<Optimization::BrentSearch>(f2, 0.0, 2.0 * PI);

    EXPECT_NEAR(x1.x, 4.712391081200089, ppx::EPS_SP * 1e2);
    EXPECT_NEAR(x2.x, 4.712391081200089, ppx::EPS_SP * 1e2);
    EXPECT_NEAR(x3.x, 4.712391081200089, ppx::EPS_SP * 1e2);

    auto f3 = [](double x)
    {
        return sin(x - 9.0 / 7.0);
    };

    x1 = fminbnd<Optimization::GoldenSearch>(f3, 1, 2.0 * PI);
    x2 = fminbnd<Optimization::QuadraticSearch>(f3, 1, 2.0 * PI);
    x3 = fminbnd<Optimization::BrentSearch>(f3, 1, 2.0 * PI);

    EXPECT_NEAR(x1.x, 5.998103062276137, ppx::EPS_SP * 1e2);
    EXPECT_NEAR(x3.x, 5.998103062276137, ppx::EPS_SP * 1e2);
}

TEST_F(NonLinOpt_TestCase, MultiDimension_OPT)
{
    auto f4 = [](const MatrixS<2, 1> &x)
    {
        return 3 * x[0] * x[0] + 2 * x[0] * x[1] + x[1] * x[1] - 4 * x[0] + 5 * x[1];
    };

    MatrixS<2, 1> result{2.250000572016352, -4.749999832427862};
    auto x1 = fminunc<Optimization::Powell>(f4, MatrixS<2, 1>{1, 1});
    EXPECT_LE(norminf((x1.x - result).eval()), 1.0e-3);

    auto f5 = [](const MatrixS<2, 1> &x)
    {
        auto sqr = x[0] * x[0] + x[1] * x[1];
        return x[0] * exp(-sqr) + sqr / 20.0;
    };
    result = {-0.669071393870035, 7.095255048137055e-07};
    x1 = fminunc<Optimization::Powell>(f5, MatrixS<2, 1>{1, 2});
    EXPECT_LE(norminf((x1.x - result).eval()), 1.0e-3);

    auto f6 = [](const MatrixS<2, 1> &x)
    {
        return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
    };

    auto d6 = [](const MatrixS<2, 1> &x)
    {
        MatrixS<2, 1> y;
        y[0] = -400 * (x[1] - x[0] * x[0]) * x[0] - 2 * (1 - x[0]);
        y[1] = 200 * (x[1] - x[0] * x[0]);
        return y;
    };

    x1 = fminunc<Optimization::Powell>(f6, MatrixS<2, 1>{4, -8});
    auto x2 = fminunc<Optimization::GradientDescent>(f6, d6, MatrixS<2, 1>{4, -8});
    auto x3 = fminunc<Optimization::ConjugateGradient>(f6, d6, MatrixS<2, 1>{4, -8});
    auto x4 = fminunc<Optimization::BGFS>(f6, d6, MatrixS<2, 1>{4, -8});
    result = {1, 1};
    EXPECT_EQ(x1.x, result);
    // EXPECT_EQ(x2.x, result);
    EXPECT_EQ(x3.x, result);
    EXPECT_EQ(x4.x, result);

    auto f7 = [](const MatrixS<2, 1> &x)
    {
        return (x[0] + 2 * x[1] - 7) * (x[0] + 2 * x[1] - 7) + (2 * x[0] + x[1] - 5) * (2 * x[0] + x[1] - 5);
    };

    auto d7 = [](const MatrixS<2, 1> &x)
    {
        MatrixS<2, 1> y;
        y[0] = 10 * x[0] + 8 * x[1] - 34;
        y[1] = 8 * x[0] + 10 * x[1] - 38;
        return y;
    };

    x1 = fminunc<Optimization::Powell>(f7, MatrixS<2, 1>{9, 8});
    x2 = fminunc<Optimization::GradientDescent>(f7, d7, MatrixS<2, 1>{9, 8});
    x3 = fminunc<Optimization::ConjugateGradient>(f7, d7, MatrixS<2, 1>{9, 8});
    x4 = fminunc<Optimization::BGFS>(f7, d7, MatrixS<2, 1>{9, 8});
    EXPECT_EQ(x1.x, x2.x);
    EXPECT_EQ(x2.x, x3.x);
    EXPECT_EQ(x3.x, x4.x);
}