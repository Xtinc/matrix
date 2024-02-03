#include "robotics.hpp"
#include "plog.h"
#include <random>
#include <fstream>
#include <chrono>

using namespace ppx;
using namespace ppx::details;

#define EXPECT_EQ(A, B) assert(A == B)
#define EXPECT_NEAR(A, B, C) assert(fabs(A - B) < C)

template <typename T>
inline void PRINT_SINGLE_ELEMENTS(const T &coll, const std::string &optcsrt = "")
{
    LOG_CH(001) << optcsrt << coll;
}

template <typename T>
inline void PRINT_LISTED_ELEMENTS(const T &coll, const std::string &optcsrt = "")
{
    LOG_CH(001) << optcsrt << coll;
}

void test_expr()
{
    MatrixS<8, 4> m;
    ones(m);
    PRINT_SINGLE_ELEMENTS(m, "m = ");
    // elem 0
    m.sub<2, 2, false>(0, 0) = {3, 3, 3};
    m.sub<3, 1, false>(3, 3) = Abs(m.sub<3, 1>(3, 3) * 2 - MatrixS<3, 1>{});
    PRINT_SINGLE_ELEMENTS(m, "m = ");
}

void test_matrix()
{
    MatrixS<4, 4> a = {1, 2, 3, 4, 5, 6, 7, 8};
    PRINT_SINGLE_ELEMENTS(a);
    PRINT_SINGLE_ELEMENTS(a(maxloc(a).first, maxloc(a).second), "max of a : ");
    a.sub<2, 2>(2, 2) = MatrixS<2, 2>{1, 1, 1, 1} + 1;
    a.sub<2, 2>(0, 0) = a.sub<2, 2>(1, 1);
    a *= -1;
    PRINT_SINGLE_ELEMENTS(a);
    PRINT_SINGLE_ELEMENTS(Abs(a).eval());
    // a(0, {0, -1}) = {89.0, 23.0, 44.0, 9.8};
    PRINT_SINGLE_ELEMENTS(a, "a = ");
    PRINT_LISTED_ELEMENTS(a, "a in list: ");
    PRINT_SINGLE_ELEMENTS(determinant(a), "determinant(a) = ");
    PRINT_SINGLE_ELEMENTS(inverse(a), "inverse(a) = ");
    MatrixS<4, 4> A(std::move(a));
    PRINT_LISTED_ELEMENTS(a, "Move Semantics, a = ");
    PRINT_SINGLE_ELEMENTS(A, "Move Semantics, A = ");
    MatrixS<2, 2> b = {1, 2, 4, 3};
    PRINT_SINGLE_ELEMENTS(inverse(b), "inverse(b) = ");
    MatrixS<1, 1> c = {3};
    PRINT_SINGLE_ELEMENTS(c.I(), "inverse(c) = ");
    MatrixS<30, 30> M{};
    std::default_random_engine eni;
    std::uniform_real_distribution<> uf(-5, 5);
    for (auto &ele : M)
    {
        ele = uf(eni);
    }
    PRINT_SINGLE_ELEMENTS(M, "M = ");
    PRINT_SINGLE_ELEMENTS(determinant(M), "determinant(M) = ");
    PRINT_SINGLE_ELEMENTS(adjugate(M), "adjoint(M) = ");
    MatrixS<4, 4> x = {1, 2, 3, 4, 5, 6, 7, 8};
    MatrixS<4, 4> y = x.T();
    MatrixS<4, 4> z;
    z.sub<4, 4>(0, 0) = x * 2 + y * 2 + 3.3 + x * y;
    PRINT_SINGLE_ELEMENTS(z, "2*(x + x^T) + 3.3 +x*x^T = ");
    MatrixS<3, 3> so3mat = {0, 3, -2, -3, 0, 1, 2, -1, 0};
    PRINT_SINGLE_ELEMENTS(vee(so3mat).exp(), "exp(s) = ");
    SO3 SO3mat = {0, 1, 0, 0, 0, 1, 1, 0, 0};
    PRINT_SINGLE_ELEMENTS(SO3mat.log(), "log(S) = ");
    PRINT_SINGLE_ELEMENTS(SE3(SO3mat, {1, 1, 1}), "T = ");
    SE3 SE3mat = {1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 3, 1};
    PRINT_SINGLE_ELEMENTS(SE3mat, "SE3mat = ");
    PRINT_SINGLE_ELEMENTS(SE3mat.log(), "log(s) = ");
    PRINT_SINGLE_ELEMENTS(hat(SE3mat.log()));
    PRINT_SINGLE_ELEMENTS(se3{1.5708, 0.0, 0.0, 0.0, 2.3562, 2.3562}.exp());
    MatrixS<4, 3> X{3.5, 1.6, 3.7, 4.3,
                    2.7, -5.7, -0.8, -9.8,
                    -3.1, -6.0, 1.9, 6.9};
    MatrixS<4, 1> Y{1, 1, 1, 1};
    PRINT_SINGLE_ELEMENTS(linsolve<Factorization::QR>(X, Y), "solved by QR = ");
    PRINT_SINGLE_ELEMENTS(linsolve<Factorization::SVD>(X, Y), "solved by SVD = ");
    PRINT_SINGLE_ELEMENTS(X, " X = ");
    MatrixS<4, 3> u{1, 4, 7, 11,
                    2, 5, 8, 1,
                    3, 6, 9, 5};
    SVD<4, 3> ru(u);
    PRINT_SINGLE_ELEMENTS(ru.u, "U = ");
    PRINT_SINGLE_ELEMENTS(ru.w, "S = ");
    PRINT_SINGLE_ELEMENTS(ru.v, "V = ");
    PRINT_SINGLE_ELEMENTS(ru.u * ru.w.diag() * ru.v.T(), "U*S*V = ");

    MatrixS<4, 4> g{2, -1, 1, 4,
                    -1, 2, -1, 1,
                    1, -1, 2, -2,
                    4, 1, -2, 3};
    PRINT_SINGLE_ELEMENTS(g, "g = ");
    EigenValue<4> eig(g);
    PRINT_SINGLE_ELEMENTS(eig.z, "eigen vector of g : ");
    PRINT_SINGLE_ELEMENTS(eig.d, "eigen value of g : ");
    // Matrix<3, 3> TA{1, 2, 3, 9, 8, 7, 5, 6, 4};
    // Unsymmeig<3> eigsolver;
    // eigsolver(TA);
    // PRINT_LISTED_ELEMENTS(eigsolver.wri, "eigvalue: ");
    // PRINT_SINGLE_ELEMENTS(eigsolver.zz, "eigen_vec: ");
}

void test_linear()
{
    std::default_random_engine eni((unsigned)time(NULL));
    std::uniform_real_distribution<> uf(-10000, 10000);
    LOG_FMT(001, "Test linear linsolver LU\n");
    for (size_t i = 0; i < 50; i++)
    {
        MatrixS<100, 100> A;
        for (auto &ele : A)
        {
            ele = uf(eni);
        }
        MatrixS<100, 1> x;
        for (auto &ele : x)
        {
            ele = uf(eni);
        }
        auto b = A * x;
        auto result = linsolve<Factorization::LU>(A, b);
        auto residual = norm2<100, 1>(result.x - x);
        PRINT_SINGLE_ELEMENTS(residual, "residual = ");
    }
    LOG_FMT(001, "Test linear linsolver SVD\n");
    for (size_t i = 0; i < 50; i++)
    {
        MatrixS<100, 100> A;
        for (auto &ele : A)
        {
            ele = uf(eni);
        }
        MatrixS<100, 1> x;
        for (auto &ele : x)
        {
            ele = uf(eni);
        }
        auto b = A * x;
        auto result = linsolve<Factorization::SVD>(A, b);
        auto residual = norm2<100, 1>(result.x - x);
        PRINT_SINGLE_ELEMENTS(residual, "residual = ");
    }
    LOG_FMT(001, "Test linear linsolver QR\n");
    for (size_t i = 0; i < 50; i++)
    {
        MatrixS<100, 100> A;
        for (auto &ele : A)
        {
            ele = uf(eni);
        }
        MatrixS<100, 1> x;
        for (auto &ele : x)
        {
            ele = uf(eni);
        }
        auto b = A * x;
        auto result = linsolve<Factorization::QR>(A, b);
        auto residual = norm2<100, 1>(result.x - x);
        PRINT_SINGLE_ELEMENTS(residual, "residual = ");
    }
}

void test_nonlinear()
{
    LOG_FMT(001, "test nonlinear 1D:\n");
    auto f1 = [](double x)
    {
        double y = 0.0;
        for (int k = -10; k < 11; k++)
        {
            y += (k + 1) * (k + 1) * cos(k * x) * exp(-1 * k * k / 2.0);
        }
        return y;
    };

    PRINT_SINGLE_ELEMENTS(fminbnd<Optimization::GoldenSearch>(f1, 1.0, 3.0), "f1 by GoldenSearch: ");
    PRINT_SINGLE_ELEMENTS(fminbnd<Optimization::QuadraticSearch>(f1, 1.0, 3.0), "f1 by QuadraticSearch: ");
    PRINT_SINGLE_ELEMENTS(fminbnd<Optimization::BrentSearch>(f1, 1.0, 3.0), "f1 by BrentSearch: ");

    auto f2 = [](double x)
    {
        return sin(x);
    };

    PRINT_SINGLE_ELEMENTS(fminbnd<Optimization::GoldenSearch>(f2, 0.0, 2.0 * PI), "f2 by GoldenSearch: ");
    PRINT_SINGLE_ELEMENTS(fminbnd<Optimization::QuadraticSearch>(f2, 0.0, 2.0 * PI), "f2 by QuadraticSearch: ");
    PRINT_SINGLE_ELEMENTS(fminbnd<Optimization::BrentSearch>(f2, 0.0, 2.0 * PI), "f2 by BrentSearch: ");

    auto f3 = [](double x)
    {
        return sin(x - 9.0 / 7.0);
    };

    PRINT_SINGLE_ELEMENTS(fminbnd<Optimization::GoldenSearch>(f3, 1, 2.0 * PI), "f3 by GoldenSearch: ");
    PRINT_SINGLE_ELEMENTS(fminbnd<Optimization::QuadraticSearch>(f3, 1, 2.0 * PI), "f3 by QuadraticSearch: ");
    PRINT_SINGLE_ELEMENTS(fminbnd<Optimization::BrentSearch>(f3, 1, 2.0 * PI), "f3 by BrentSearch: ");

    LOG_FMT(001, "test nonlinear ND:\n");
    auto f4 = [](const MatrixS<2, 1> &x)
    {
        return 3 * x[0] * x[0] + 2 * x[0] * x[1] + x[1] * x[1] - 4 * x[0] + 5 * x[1];
    };

    PRINT_SINGLE_ELEMENTS(fminunc<Optimization::Powell>(f4, MatrixS<2, 1>{1, 1}), "f4 by Powell: ");

    auto f5 = [](const MatrixS<2, 1> &x)
    {
        auto sqr = x[0] * x[0] + x[1] * x[1];
        return x[0] * x[1] * exp(-sqr) + sqr / 20.0;
    };

    PRINT_SINGLE_ELEMENTS(fminunc<Optimization::Powell>(f5, MatrixS<2, 1>{1, 2}), "f5 by Powell: ");

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

    PRINT_SINGLE_ELEMENTS(fminunc<Optimization::Powell>(f6, MatrixS<2, 1>{4, -8}), "f6 by Powell: ");
    PRINT_SINGLE_ELEMENTS(fminunc<Optimization::GradientDescent>(f6, d6, MatrixS<2, 1>{4, -8}), "f6 by GradientDescent: ");
    PRINT_SINGLE_ELEMENTS(fminunc<Optimization::ConjugateGradient>(f6, d6, MatrixS<2, 1>{4, -8}), "f6 by ConjuateGradient: ");
    PRINT_SINGLE_ELEMENTS(fminunc<Optimization::BGFS>(f6, d6, MatrixS<2, 1>{4, -8}), "f6 by BGFS: ");

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

    PRINT_SINGLE_ELEMENTS(fminunc<Optimization::Powell>(f7, MatrixS<2, 1>{9, 8}), "f7 by Powell: ");
    PRINT_SINGLE_ELEMENTS(fminunc<Optimization::GradientDescent>(f7, d7, MatrixS<2, 1>{9, 8}), "f7 by GradientDescent: ");
    PRINT_SINGLE_ELEMENTS(fminunc<Optimization::ConjugateGradient>(f7, d7, MatrixS<2, 1>{9, 8}), "f7 by ConjuateGradient: ");
    PRINT_SINGLE_ELEMENTS(fminunc<Optimization::BGFS>(f7, d7, MatrixS<2, 1>{9, 8}), "f7 by BGFS: ");
}

void test_statics()
{
    MatrixS<2, 1> mu{-2, 8};
    MatrixS<2, 2> sigma{0.4, 0.0, 0.0, 1.2};
    MatrixS<2, 1> x;
    MultiGaussianDistribution<2> d(mu, sigma);
    PRINT_SINGLE_ELEMENTS(d.pdf(x), "pdf of [0,0] = ");
    PRINT_SINGLE_ELEMENTS(d.pdf(MatrixS<2, 1>{-0.6, -0.6}), "pdf of [-0.6,-0.6] = ");
    std::vector<MatrixS<2, 1>> x1;
    MultiGaussianDistribution<2> d2(MatrixS<2, 1>{2, 1}, MatrixS<2, 2>{4, 1, 1, 4});
    MultiGaussianDistribution<2> d3(MatrixS<2, 1>{6, 4}, MatrixS<2, 2>{0.25, 1.6, 2.0, 16.0});
    MixedGaussianDistribution<2, 3> mg;
    mg.setcomp(0, d, 0.3);
    mg.setcomp(1, d2, 0.5);
    mg.setcomp(2, d3, 0.2);
    constexpr int ITMAX = 1000;
    for (int n = 0; n < ITMAX; ++n)
    {
        x1.emplace_back(mg());
    }
    PRINT_SINGLE_ELEMENTS(avg(x1), "mean = ");
    // for (size_t i = 0; i < ITMAX; i++)
    // {
    //     PRINT_LISTED_ELEMENTS(x1[i]);
    // }

    MixedGaussianDistribution<2, 3> mge;
    mge.setcomp(0, MultiGaussianDistribution<2>({1, 1}, {1, 0, 0, 1}), 0.5);
    mge.setcomp(1, MultiGaussianDistribution<2>({2, 2}, {1, 0, 0, 1}), 0.2);
    mge.setcomp(2, MultiGaussianDistribution<2>({4, 4}, {1, 0, 0, 1}), 0.3);
    mge.fit(x1);
    PRINT_SINGLE_ELEMENTS(mge);
    // auto e1 = mean(x1);
    // auto e2 = mean(x2);
    // auto s1 = var(x1);
    // auto s2 = var(x2);
    // PRINT_SINGLE_ELEMENTS(e1, " mean of x1 = ");
    // PRINT_SINGLE_ELEMENTS(s1, " var of x1 = ");
    // PRINT_SINGLE_ELEMENTS(e2, " mean of x2 = ");
    // PRINT_SINGLE_ELEMENTS(s2, " var of x2 = ");
}

void test_lieGroup()
{
    {
        MatrixS<3, 1> vec{1, 2, 3};
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
        MatrixS<4, 4> expected{0.0, 0.0, 0.0, 0.0,
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
        MatrixS<6, 6> expected{1, 0, 0, 0, 3, 0,
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
        MatrixS<4, 4> se3mat = {0.0, 0.0, 0.0, 0.0,
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
    auto s = StatusCode::NORMAL;
    kinematics<6>::Q q;
    while (s == StatusCode::NORMAL)
    {
        random(q, -PI + 1e-3, PI - 1e-3);
        auto J = UR5.jacobiSpace(q);
        LU<6> lu(J);
        s = lu.s;
    }
    q = {0.336058, 1.69911, -0.445658, -1.25227, 0.00013744, -2.71554};
    SE3 TargetPose = UR5.forwardSpace(q);
    auto t1 = std::chrono::system_clock::now();
    for (size_t i = 0; i < 2000; i++)
    {
        q = UR5.inverseSpace(TargetPose, q - 0.05);
    }
    auto t2 = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 2000.0 << " us per iter has elapsed.\n";
    PRINT_SINGLE_ELEMENTS(q, "IKSpace = ");
    std::cout << "singualr q:" << q << "\n";
}

int main(int, char **)
{
    ppx::initialize_log(opts);

    test_expr();
    test_matrix();
    test_linear();
    test_statics();
    test_lieGroup();
    test_nonlinear();
    test_robotics();
    return EXIT_SUCCESS;
}
