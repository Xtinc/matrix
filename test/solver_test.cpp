#include "gtest/gtest.h"
#include "liegroup.hpp"

using namespace ppx;

class LinEqn_TestCase : public ::testing::Test
{
public:
    LinEqn_TestCase() : eni(std::random_device()()){};
    std::default_random_engine eni;
};

TEST_F(LinEqn_TestCase, Quadtratic)
{
    auto res = quadsolve(1.0, -2 * 12345678, -1);
    EXPECT_EQ(res.s, StatusCode::NORMAL);
    EXPECT_NEAR(res.x[1] / -4.050000332100021e-08, 1, 1e-7);
}

TEST_F(LinEqn_TestCase, MGS)
{
    MatrixS<3, 3> A{
        {{1.0001, 0.5000, 0.3333},
         {0.5000, 0.3334, 0.2500},
         {0.3333, 0.2500, 0.2001}},
        Ori::Row};
    auto Q = MGS(A);
    auto residual = norm1<3, 3>(Q.T() * Q - eye<3>());
    EXPECT_TRUE(residual < EPS_SP);
    MatrixS<100, 100> B;
    for (size_t i = 0; i < 100; i++)
    {
        random(B);
        auto R = MGS(B);
        residual = norm1<3, 3>(R.T() * R - eye<100>());
        EXPECT_TRUE(residual < EPS_SP);
    }
}

TEST_F(LinEqn_TestCase, LU_decompose)
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
        auto result = linsolve<Factorization::LU>(A, b);
        if (result.s != StatusCode::SINGULAR)
        {
            auto residual = norm2(result.x - x);
            EXPECT_LE(residual, ppx::EPS_SP * 1e3);
        }
    }
}

TEST_F(LinEqn_TestCase, QR_decompose)
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
        auto result = linsolve<Factorization::QR>(A, b);
        if (result.s != StatusCode::SINGULAR)
        {
            auto residual = norm2(result.x - x);
            EXPECT_LE(residual, ppx::EPS_SP * 1e3);
        }
    }
}

TEST_F(LinEqn_TestCase, Matrix_Norm2)
{
    MatrixS<5, 4> A{{0.814723686393179,
                     0.905791937075619,
                     0.126986816293506,
                     0.913375856139019,
                     0.632359246225410},
                    {0.0975404049994095,
                     0.278498218867048,
                     0.546881519204984,
                     0.957506835434298,
                     0.964888535199277},
                    {0.157613081677548,
                     0.970592781760616,
                     0.957166948242946,
                     0.485375648722841,
                     0.800280468888800},
                    {0.141886338627215,
                     0.421761282626275,
                     0.915735525189067,
                     0.792207329559554,
                     0.959492426392903}};
    EXPECT_NEAR(norm2(A), 2.993584186569183, 10 * EPS_DP);
    MatrixS<4, 5> B{
        {0.655740699156587,
         0.0357116785741896,
         0.849129305868777,
         0.933993247757551},
        {0.678735154857774,
         0.757740130578333,
         0.743132468124916,
         0.392227019534168},
        {0.655477890177557,
         0.171186687811562,
         0.706046088019609,
         0.0318328463774207},
        {0.276922984960890,
         0.0461713906311539,
         0.0971317812358475,
         0.823457828327293},
        {0.694828622975817,
         0.317099480060861,
         0.950222048838355,
         0.0344460805029088}};
    EXPECT_NEAR(norm2(B), 2.386200992922722, 10 * EPS_DP);
    MatrixS<4, 4> C{
        {0.438744359656398,
         0.381558457093008,
         0.765516788149002,
         0.795199901137063},
        {0.186872604554379,
         0.489764395788231,
         0.445586200710900,
         0.646313010111265},
        {0.709364830858073,
         0.754686681982361,
         0.276025076998578,
         0.679702676853675},
        {0.655098003973841,
         0.162611735194631,
         0.118997681558377,
         0.498364051982143}};
    EXPECT_NEAR(norm2(C), 2.075457904661895, 10 * EPS_DP);
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
            auto residual = norm2(result.x - x);
            EXPECT_LE(residual, ppx::EPS_SP * 1e3);
        }
    }
}