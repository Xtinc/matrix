#include "gtest/gtest.h"
#include "liegroup.hpp"

using namespace ppx;

class MatrixS_TestCase : public ::testing::Test
{
public:
    MatrixS_TestCase() = default;
};

TEST_F(MatrixS_TestCase, ctor)
{
    MatrixS<1, 1> A = {};
    EXPECT_EQ(A[0], 0);
    MatrixS<2, 2> B{1, 2, 3, 4};
    EXPECT_EQ(B(1, 1), 4);
    MatrixS<2, 2> C(std::array<int, 4>{1, 2, 3, 4});
    MatrixS<2, 2> D(std::vector<int>{1, 2, 3, 4});
    EXPECT_EQ(B, C);
    EXPECT_EQ(C, D);
    EXPECT_EQ(sizeof(MatrixS<2, 2>), (sizeof(double)) * 4);
    EXPECT_EQ(sizeof(MatrixS<20, 20>), sizeof(std::vector<double>));
}

TEST_F(MatrixS_TestCase, expr)
{
    MatrixS<4, 4> a = {1, 2, 3, 4, 5, 6, 7, 8};
    a.sub<4, 2, false>(0, 2) = {1, 2, 3, 4, 5, 6, 7, 8};
    MatrixS<4, 4> result{1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8};
    EXPECT_EQ(a, result);
    a.sub<4, 1, false>(0, 0) = a.sub<4, 1>(0, 0) - a.sub<4, 1>(0, 2);
    a.sub<4, 1, false>(0, 1) = a.sub<4, 1>(0, 1) - a.sub<4, 1>(0, 3);
    a.sub<4, 2>(0, 2) = a.sub<4, 2>(0, 2) - a.sub<4, 2>(0, 2);
    a.sub<4, 2>(0, 2) = a.sub<4, 2>(0, 2) - a.sub<4, 2>(0, 2);
    EXPECT_EQ(a, decltype(a)::zeros());
    MatrixS<3, 3> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    MatrixS<3, 3> y = x.T();
    MatrixS<3, 3> z;
    z.sub<3, 3>(0, 0) = (x * 2 + y * 2 + x * y) / 2 - MatrixS<3, 3>{35.0, 45.0, 55.0, 45.0, 56.5, 68.0, 55.0, 68.0, 81.0};
    EXPECT_EQ(z, decltype(z)::zeros());
}

TEST_F(MatrixS_TestCase, func)
{
    MatrixS<2, 2> x = {1, 2, 4, 3};
    MatrixS<2, 2> result{-0.6, 0.4, 0.8, -0.2};
    EXPECT_NEAR(x.det(), -5.0, ppx::EPS_DP);
    EXPECT_EQ(x.I(), result);
    MatrixS<1, 1> c = {3};
    EXPECT_NEAR(c.I()[0], 1.0 / 3.0, ppx::EPS_DP);
}

TEST_F(MatrixS_TestCase, fac)
{
    MatrixS<4, 3> X{3.5, 1.6, 3.7, 4.3,
                    2.7, -5.7, -0.8, -9.8,
                    -3.1, -6.0, 1.9, 6.9};
    MatrixS<4, 1> Y{1, 1, 1, 1};
    MatrixS<3, 1> result{0.271349846985587, -0.030388140654028, -0.062228084118990};
    auto x1 = linsolve<Factorization::QR>(X, Y);
    auto x2 = linsolve<Factorization::SVD>(X, Y);
    EXPECT_EQ(x1.x, result);
    EXPECT_EQ(x2.x, result);

    MatrixS<4, 3> u{1, 4, 7, 11,
                    2, 5, 8, 1,
                    3, 6, 9, 5};
    MatrixS<3, 1> v{};
    MatrixS<3, 3> w{};
    auto ru = svdcmp(u, v, w);
    EXPECT_EQ(ru.x * v.diag() * w.T(), u);

    MatrixS<4, 4> g{2, -1, 1, 4,
                    -1, 2, -1, 1,
                    1, -1, 2, -2,
                    4, 1, -2, 3};
    auto e1 = eig<EigenSystem::SymValAndVecSorted>(g);
    MatrixS<4, 1> result2{-2.674416988218162, 1.0, 3.945104914279287, 6.729312073938871};
    EXPECT_EQ(e1.val, result2);
    MatrixS<4, 4> result3{-0.657702449551592, 0.127000127000191, 0.452113001308500, 0.588975627376541,
                          -0.196417934159944, 0.762000762001143, -0.611122207034796, 0.0854662618752970,
                          0.367242010902251, 0.635000635000952, 0.642924048180056, -0.220354639725665,
                          0.627678889578617, 0.0, -0.0936597586400323, 0.772817611852141};
    result3 = result3.T() * (-1);
    EXPECT_EQ(e1.vec, result3);
}

TEST_F(MatrixS_TestCase, cat)
{
    MatrixS<2, 2> x = {1, 2, 3, 4};
    MatrixS<2, 2> y = {5, 6, 7, 8};
    MatrixS<2, 4> expect1{1, 2, 3, 4, 5, 6, 7, 8};
    MatrixS<4, 2> expect2{1, 2, 5, 6, 3, 4, 7, 8};
    EXPECT_EQ(expect1, concat<Orientation::Horizontal>(x, y));
    EXPECT_EQ(expect2, concat<Orientation::Vertical>(x, y));

    MatrixS<2, 3> a{1, 2, 3, 4, 5, 6};
    MatrixS<2, 1> b{9, 9};
    MatrixS<1, 3> c{-1, -1, -1};
    MatrixS<2, 4> expect3{1, 2, 3, 4, 5, 6, 9, 9};
    MatrixS<3, 3> expect4{1, 2, -1, 3, 4, -1, 5, 6, -1};
    EXPECT_EQ(expect3, concat(a, b));
    EXPECT_EQ(expect4, concat(a, c));
    MatrixS<3, 1> expect5{9, 9, 1};
    MatrixS<3, 1> expect6{1, 9, 9};
    MatrixS<1, 4> expect7{-1, -1, -1, 3};
    EXPECT_EQ(expect5, concat(b, 1));
    EXPECT_EQ(expect6, concat(1, b));
    EXPECT_EQ(expect7, concat(c, 3));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}