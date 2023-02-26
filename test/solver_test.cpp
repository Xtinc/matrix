#include "gtest/gtest.h"
#include "liegroup.hpp"

using namespace ppx;

class LinEqn_TestCase : public ::testing::Test
{
public:
    LinEqn_TestCase() : eni(std::random_device()()){};
    std::default_random_engine eni;
};

TEST_F(LinEqn_TestCase, LU_decompose)
{
    std::uniform_real_distribution<> uf(-1.0e5, 1.0e5);
    for (size_t i = 0; i < 1000; i++)
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
        if (result.s != StatusCode::SINGULAR)
        {
            auto residual = norm2((result.x - x).eval());
            EXPECT_LE(residual, ppx::EPS_SP * 1e3);
        }
    }
}

TEST_F(LinEqn_TestCase, QR_decompose)
{
    std::uniform_real_distribution<> uf(-1.0e5, 1.0e5);
    for (size_t i = 0; i < 1000; i++)
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
        if (result.s != StatusCode::SINGULAR)
        {
            auto residual = norm2((result.x - x).eval());
            EXPECT_LE(residual, ppx::EPS_SP * 1e3);
        }
    }
}

TEST_F(LinEqn_TestCase, SVD_decompose)
{
    std::uniform_real_distribution<> uf(-1.0e5, 1.0e5);
    for (size_t i = 0; i < 100; i++)
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
        if (result.s != StatusCode::SINGULAR)
        {
            auto residual = norm2((result.x - x).eval());
            EXPECT_LE(residual / norm2(x), 1);
        }
    }
}