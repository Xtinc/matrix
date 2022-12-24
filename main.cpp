#include "robotics.hpp"
#include "ppxlog.h"
#include <random>
#include <cassert>

using namespace ppx;
using namespace ppx::details;

#define EXPECT_EQ(A, B) assert(A == B)
#define EXPECT_NEAR(A, B, C) assert(fabs(A - B) < C)

void test_matrix()
{
    Matrix<4, 4> a = {1, 2, 3, 4, 5, 6, 7, 8};
    PRINT_SINGLE_ELEMENTS(a);
    a({2, 3}, {2, 3}) = Matrix<2, 2>{1, 1, 1, 1};
    a(0, {0, -1}) = {89.0, 23.0, 44.0, 9.8};
    PRINT_SINGLE_ELEMENTS(a, "a = ");
    PRINT_LISTED_ELEMENTS(a, "a in list: ");
    PRINT_SINGLE_ELEMENTS(determinant(a), "determinant(a) = ");
    PRINT_SINGLE_ELEMENTS(inverse(a), "inverse(a) = ");
    Matrix<4, 4> A(std::move(a));
    PRINT_LISTED_ELEMENTS(a, "Move Sementics, a = ");
    PRINT_SINGLE_ELEMENTS(A, "Move Sementics, A = ");
    Matrix<2, 2> b = {1, 2, 4, 3};
    PRINT_SINGLE_ELEMENTS(inverse(b), "inverse(b) = ");
    Matrix<1, 1> c = {3};
    PRINT_SINGLE_ELEMENTS(c.I(), "inverse(c) = ");
    Matrix<30, 30> M{};
    std::default_random_engine eni;
    std::uniform_real_distribution<> uf(-5, 5);
    for (auto &ele : M)
    {
        ele = uf(eni);
    }
    PRINT_SINGLE_ELEMENTS(M, "M = ");
    PRINT_SINGLE_ELEMENTS(determinant(M), "determinant(M) = ");
    PRINT_SINGLE_ELEMENTS(adjugate(M), "adjoint(M) = ");
    Matrix<4, 4> x = {1, 2, 3, 4, 5, 6, 7, 8};
    Matrix<4, 4> y = x.T();
    Matrix<4, 4> z = x * 2 + y * 2 + 3.3 + x * y;
    PRINT_SINGLE_ELEMENTS(z, "2*(x + x^T) + 3.3 +x*x^T = ");
    Matrix<3, 3> so3mat = {0, 3, -2, -3, 0, 1, 2, -1, 0};
    PRINT_SINGLE_ELEMENTS(vee(so3mat).exp(), "exp(s) = ");
    SO3 SO3mat = {0, 1, 0, 0, 0, 1, 1, 0, 0};
    PRINT_SINGLE_ELEMENTS(SO3mat.log(), "log(S) = ");
    PRINT_SINGLE_ELEMENTS(SE3(SO3mat, {1, 1, 1}), "T = ");
    SE3 SE3mat = {1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 3, 1};
    PRINT_SINGLE_ELEMENTS(SE3mat, "SE3mat = ");
    PRINT_SINGLE_ELEMENTS(SE3mat.log(), "log(s) = ");
    PRINT_SINGLE_ELEMENTS(hat(SE3mat.log()));
    PRINT_SINGLE_ELEMENTS(se3{1.5708, 0.0, 0.0, 0.0, 2.3562, 2.3562}.exp());
    Matrix<4, 3> X{3.5, 1.6, 3.7, 4.3,
                   2.7, -5.7, -0.8, -9.8,
                   -3.1, -6.0, 1.9, 6.9};
    Matrix<4, 1> Y{1, 1, 1, 1};
    bool sing = true;
    PRINT_SINGLE_ELEMENTS(solve<factorization::QR>(X, Y, sing), "QR    x = ");
    PRINT_SINGLE_ELEMENTS(solve<factorization::SVD>(X, Y, sing), "SVD  x = ");
    PRINT_SINGLE_ELEMENTS(X, " X = ");
    Matrix<4, 3> u{1, 4, 7, 11,
                   2, 5, 8, 1,
                   3, 6, 9, 5};
    Matrix<3, 1> v{};
    Matrix<3, 3> w{};
    u = svdcmp(u, v, w, sing);
    PRINT_SINGLE_ELEMENTS(u, "U = ");
    PRINT_SINGLE_ELEMENTS(v, "S = ");
    PRINT_SINGLE_ELEMENTS(w, "V = ");
    PRINT_SINGLE_ELEMENTS(u * Matrix<3, 3>::diag({v[0], v[1], v[2]}) * w.T(), "U*S*V = ");
}

void test_linear()
{
    std::default_random_engine eni((unsigned)time(NULL));
    std::uniform_real_distribution<> uf(-10000, 10000);
    printf("Test linear solver LU\n");
    for (size_t i = 0; i < 50; i++)
    {
        Matrix<100, 100> A;
        for (auto &ele : A)
        {
            ele = uf(eni);
        }
        Matrix<100, 1> x;
        for (auto &ele : x)
        {
            ele = uf(eni);
        }
        auto b = A * x;
        bool sing = false;
        auto result = solve<factorization::LU>(A, b, sing);
        auto residual = norm2(Matrix<100, 1>(result - x));
        PRINT_SINGLE_ELEMENTS(residual, "residual = ");
    }
    printf("Test linear solver SVD\n");
    for (size_t i = 0; i < 50; i++)
    {
        Matrix<100, 100> A;
        for (auto &ele : A)
        {
            ele = uf(eni);
        }
        Matrix<100, 1> x;
        for (auto &ele : x)
        {
            ele = uf(eni);
        }
        auto b = A * x;
        bool sing = false;
        auto result = solve<factorization::SVD>(A, b, sing);
        auto residual = norm2(Matrix<100, 1>(result - x));
        PRINT_SINGLE_ELEMENTS(residual, "residual = ");
    }
    printf("Test linear solver QR\n");
    for (size_t i = 0; i < 50; i++)
    {
        Matrix<100, 100> A;
        for (auto &ele : A)
        {
            ele = uf(eni);
        }
        Matrix<100, 1> x;
        for (auto &ele : x)
        {
            ele = uf(eni);
        }
        auto b = A * x;
        bool sing = false;
        auto result = solve<factorization::QR>(A, b, sing);
        auto residual = norm2(Matrix<100, 1>(result - x));
        PRINT_SINGLE_ELEMENTS(residual, "residual = ");
    }
}

void test_nonlinear()
{
    printf("test nonlinear 1D:\n");
    auto f1 = [](double x)
    {
        double y = 0.0;
        for (int k = -10; k < 11; k++)
        {
            y += (k + 1) * (k + 1) * cos(k * x) * exp(-1 * k * k / 2.0);
        }
        return y;
    };

    PRINT_SINGLE_ELEMENTS(fminbnd<optimization::GoldenSearch>(f1, 1.0, 3.0), "f1 by GoldenSearch: ");
    PRINT_SINGLE_ELEMENTS(fminbnd<optimization::QuadraticSearch>(f1, 1.0, 3.0), "f1 by QuadraticSearch: ");
    PRINT_SINGLE_ELEMENTS(fminbnd<optimization::BrentSearch>(f1, 1.0, 3.0), "f1 by BrentSearch: ");

    auto f2 = [](double x)
    {
        return sin(x);
    };

    PRINT_SINGLE_ELEMENTS(fminbnd<optimization::GoldenSearch>(f2, 0.0, 2.0 * gl_rep_pi), "f2 by GoldenSearch: ");
    PRINT_SINGLE_ELEMENTS(fminbnd<optimization::QuadraticSearch>(f2, 0.0, 2.0 * gl_rep_pi), "f2 by QuadraticSearch: ");
    PRINT_SINGLE_ELEMENTS(fminbnd<optimization::BrentSearch>(f2, 0.0, 2.0 * gl_rep_pi), "f2 by BrentSearch: ");

    auto f3 = [](double x)
    {
        return sin(x - 9.0 / 7.0);
    };

    PRINT_SINGLE_ELEMENTS(fminbnd<optimization::GoldenSearch>(f3, 1, 2.0 * gl_rep_pi), "f3 by GoldenSearch: ");
    PRINT_SINGLE_ELEMENTS(fminbnd<optimization::QuadraticSearch>(f3, 1, 2.0 * gl_rep_pi), "f3 by QuadraticSearch: ");
    PRINT_SINGLE_ELEMENTS(fminbnd<optimization::BrentSearch>(f3, 1, 2.0 * gl_rep_pi), "f3 by BrentSearch: ");

    printf("test nonlinear ND:\n");
    auto f4 = [](const Matrix<2, 1> &x)
    {
        return 3 * x[0] * x[0] + 2 * x[0] * x[1] + x[1] * x[1] - 4 * x[0] + 5 * x[1];
    };

    PRINT_SINGLE_ELEMENTS(fminunc<optimization::Powell>(f4, Matrix<2, 1>{1, 1}), "f4 by Powell: ");

    auto f5 = [](const Matrix<2, 1> &x)
    {
        auto sqr = x[0] * x[0] + x[1] * x[1];
        return x[0] * x[1] * exp(-sqr) + sqr / 20.0;
    };

    PRINT_SINGLE_ELEMENTS(fminunc<optimization::Powell>(f5, Matrix<2, 1>{1, 2}), "f5 by Powell: ");

    auto f6 = [](const Matrix<2, 1> &x)
    {
        return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
    };

    auto d6 = [](const Matrix<2, 1> &x)
    {
        Matrix<2, 1> y;
        y[0] = -400 * (x[1] - x[0] * x[0]) * x[0] - 2 * (1 - x[0]);
        y[1] = 200 * (x[1] - x[0] * x[0]);
        return y;
    };

    PRINT_SINGLE_ELEMENTS(fminunc<optimization::Powell>(f6, Matrix<2, 1>{4, -8}), "f6 by Powell: ");
    PRINT_SINGLE_ELEMENTS(fminunc<optimization::GradientDescent>(f6, d6, Matrix<2, 1>{4, -8}), "f6 by GradientDescent: ");
    PRINT_SINGLE_ELEMENTS(fminunc<optimization::ConjuateGradient>(f6, d6, Matrix<2, 1>{4, -8}), "f6 by ConjuateGradient: ");
    PRINT_SINGLE_ELEMENTS(fminunc<optimization::BGFS>(f6, d6, Matrix<2, 1>{4, -8}), "f6 by BGFS: ");

    auto f7 = [](const Matrix<2, 1> &x)
    {
        return (x[0] + 2 * x[1] - 7) * (x[0] + 2 * x[1] - 7) + (2 * x[0] + x[1] - 5) * (2 * x[0] + x[1] - 5);
    };

    auto d7 = [](const Matrix<2, 1> &x)
    {
        Matrix<2, 1> y;
        y[0] = 10 * x[0] + 8 * x[1] - 34;
        y[1] = 8 * x[0] + 10 * x[1] - 38;
        return y;
    };

    PRINT_SINGLE_ELEMENTS(fminunc<optimization::Powell>(f7, Matrix<2, 1>{9, 8}), "f7 by Powell: ");
    PRINT_SINGLE_ELEMENTS(fminunc<optimization::GradientDescent>(f7, d7, Matrix<2, 1>{9, 8}), "f7 by GradientDescent: ");
    PRINT_SINGLE_ELEMENTS(fminunc<optimization::ConjuateGradient>(f7, d7, Matrix<2, 1>{9, 8}), "f7 by ConjuateGradient: ");
    PRINT_SINGLE_ELEMENTS(fminunc<optimization::BGFS>(f7, d7, Matrix<2, 1>{9, 8}), "f7 by BGFS: ");
}

void test_lieGroup()
{
    {
        Matrix<3, 1> vec{1, 2, 3};
        SO3 expected{0, 3, -2, -3, 0, 1, 2, -1, 0};
        SO3 result = hat(vec);
        EXPECT_EQ(result, expected);
        PRINT_SINGLE_ELEMENTS(result, "exp(so3) = ");
    }
    {
        SE3 Tinput{1, 0, 0, 0,
                   0, 0, 1, 0,
                   0, -1, 0, 0,
                   0, 0, 3, 1};
        Matrix<4, 4> expected{0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 1.57079633, 0.0,
                              0.0, -1.57079633, 0.0, 0.0,
                              0.0, 2.35619449, 2.35619449, 0.0};
        auto result = hat(Tinput.log());
        EXPECT_EQ(result, expected);
        PRINT_SINGLE_ELEMENTS(result, "log(SE3) = ");
    }
    {
        SE3 T{1, 0, 0, 0,
              0, 0, 1, 0,
              0, -1, 0, 0,
              0, 0, 3, 1};
        Matrix<6, 6> expected{1, 0, 0, 0, 3, 0,
                              0, 0, 1, 0, 0, 0,
                              0, -1, 0, 3, 0, 0,
                              0, 0, 0, 1, 0, 0,
                              0, 0, 0, 0, 0, 1,
                              0, 0, 0, 0, -1, 0};
        auto result = T.Adt();
        EXPECT_EQ(expected, result);
        PRINT_SINGLE_ELEMENTS(result, "SE3.Adjoint = ");
    }
    {
        Matrix<4, 4> se3mat = {0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 1.5708, 0.0,
                               0.0, -1.5708, 0.0, 0.0,
                               0.0, 2.3562, 2.3562, 0.0};
        SE3 result{{1, 0, 0, 0, 0, 1, 0, -1, 0, 0}, {0, 0, 3}};
        auto cal = vee(se3mat).exp();
        PRINT_SINGLE_ELEMENTS(cal, "exp(se3) = ");
        for (size_t i = 0; i < 16; i++)
        {
            EXPECT_NEAR(cal[i], result[i], 1.0e-5);
        }
    }
}

void test_robotics()
{
    {
        kinematics<6> UR5;
        SE3 F6{-1.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0,
               0.0, 1.0, 0.0, 0.0,
               0.817, 0.191, -0.006, 1.0};
        UR5.setJoint<0>({"R1", se3{0, 0, 1, 0, 0, 0, 0}, SE3()});
        UR5.setJoint<1>({"R2", se3{0.0, 1.0, 0.0, -0.089, 0.0, 0.0}, SE3{}});
        UR5.setJoint<2>({"R3", se3{0.0, 1.0, 0.0, -0.089, 0.0, 0.425}, SE3{}});
        UR5.setJoint<3>({"R4", se3{0.0, 1.0, 0.0, -0.089, 0.0, 0.817}, SE3{}});
        UR5.setJoint<4>({"R5", se3{0.0, 0.0, -1.0, -0.109, 0.817, 0.0}, SE3{}});
        UR5.setJoint<5>({"R6", se3{0.0, 1.0, 0.0, 0.006, 0.0, 0.817}, F6});
        PRINT_SINGLE_ELEMENTS(UR5.forwardSpace("R6", {0.0, -0.5 * gl_rep_pi, 0.0, 0.0, 0.5 * gl_rep_pi, 0.0}), "Forward(R6) = ");
        PRINT_SINGLE_ELEMENTS(UR5.jacobiSpace({0.0, -0.5 * gl_rep_pi, 0.0, 0.0, 0.5 * gl_rep_pi, 0.0}), "Jacobi = ");
        PRINT_SINGLE_ELEMENTS(UR5.jacobiSpace(std::array<std::string, 3>{"R1", "R2", "R3"}, {0.0, -0.5 * gl_rep_pi, 0.0, 0.0, 0.5 * gl_rep_pi, 0.0}), "Jacobi(3) = ");
        SE3 TargetPose{0.0, 1.0, 0.0, 0.0,
                       -1.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 1.0, 0.0,
                       0.095, 0.109, 0.988, 1.0};
        PRINT_SINGLE_ELEMENTS(UR5.inverseSpace(TargetPose, {0.0, -1.5, 0.0, 0.0, 1.5, 0.0}), "IKSpace = ");
    }
}

int main(int, char **)
{
    ppx::initialize_log("./", "test", 10);
    LOG_INFO << "test_matrix";
    test_matrix();
    LOG_INFO << "test_linear";
    test_linear();
    LOG_INFO << "test_lieGroup";
    test_lieGroup();
    LOG_INFO << "test_robotics";
    test_robotics();
    LOG_INFO << "test_nonlinear";
    test_nonlinear();
    LOG_INFO << "test end";
}
