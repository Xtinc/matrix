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

TEST_F(NonLinOpt_TestCase, RidgeRegress)
{
    MatrixS<16, 6> D{
        1300.0, 1300.0, 1300.0, 1300.0, 1300.0, 1300.0, 1200.0, 1200.0, 1200.0, 1200.0, 1200.0, 1200.0, 1100.0, 1100.0, 1100.0, 1100.0,
        7.5, 9.0, 11.0, 13.5, 17.0, 23.0, 5.3, 7.5, 11.0, 13.5, 17.0, 23.0, 5.3, 7.5, 11.0, 17.0,
        0.012, 0.012, 0.0115, 0.013, 0.0135, 0.012, 0.040, 0.038, 0.032, 0.026, 0.034, 0.041, 0.084, 0.098, 0.092, 0.086,
        9750.0, 11700.0, 14300.0, 17550.0, 22100.0, 29900.0, 6360.0, 9000.0, 13200.0, 16200.0, 20400.0, 27600.0, 5830.0, 8250.0, 12100.0, 18700.0,
        15.6, 15.6, 14.95, 16.9, 17.55, 15.6, 48.0, 45.6, 38.4, 31.2, 40.8, 49.2, 92.4, 107.8, 101.2, 94.6,
        0.090, 0.108, 0.1265, 0.1755, 0.2295, 0.276, 0.212, 0.285, 0.352, 0.351, 0.578, 0.943, 0.4452, 0.735, 1.012, 1.462};
    MatrixS<16, 1> y{49.0, 50.2, 50.5, 48.5, 47.5, 44.5, 28.0, 31.5,
                     34.5, 35.0, 38.0, 38.5, 15.0, 17.0, 20.5, 29.5};
    for (size_t i = 0; i < 501; i++)
    {
        auto x = regress(y, D, 1e-6 * i).x;
        std::cout << x << std::endl;
    }
}