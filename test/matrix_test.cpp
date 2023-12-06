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
    MatrixS<3, 2> E{1, 2, 3, 4, 5, 6};
    MatrixS<3, 2> F{{1, 2, 3},
                    {4, 5, 6}};
    MatrixS<3, 2> G{{{1, 4}, {2, 5}, {3, 6}}, Ori::Row};
    EXPECT_EQ(E, F);
    EXPECT_EQ(F, G);
    MatrixS<3, 2> H{{1, 2, 3, 4},
                    {4, 5, 6}};
    EXPECT_EQ(G, H);
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
    EXPECT_EQ(a, decltype(a){});
    MatrixS<3, 3> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    MatrixS<3, 3> y = x.T();
    MatrixS<3, 3> z;
    z.sub<3, 3>(0, 0) =
        (x * 2 + y * 2 + x * y) / 2 - MatrixS<3, 3>{35.0, 45.0, 55.0, 45.0, 56.5, 68.0, 55.0, 68.0, 81.0};
    EXPECT_EQ(z, decltype(z){});
    MatrixS<3, 3> w = Sqrt(pwmul(x, x));
    EXPECT_EQ(w, x);
    MatrixS<3, 3> v = Abs(pwdiv(pwmul(x, x.sub<3, 3>(0, 0)), 1 + x.sub<3, 3>(0, 0) - 1));
    EXPECT_EQ(w, v);
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
    SVD<4, 3> ru(u);
    EXPECT_EQ(ru.u * ru.w.diag() * ru.v.T(), u);

    MatrixS<4, 4> g{2, -1, 1, 4,
                    -1, 2, -1, 1,
                    1, -1, 2, -2,
                    4, 1, -2, 3};
    EigenValue<4> e1(g, true);
    MatrixS<4, 1> result2{-2.674416988218162, 1.0, 3.945104914279287, 6.729312073938871};
    EXPECT_EQ(e1.d, result2);
    MatrixS<4, 4> result3{-0.657702449551592, 0.127000127000191, 0.452113001308500, 0.588975627376541,
                          -0.196417934159944, 0.762000762001143, -0.611122207034796, 0.0854662618752970,
                          0.367242010902251, 0.635000635000952, 0.642924048180056, -0.220354639725665,
                          0.627678889578617, 0.0, -0.0936597586400323, 0.772817611852141};
    result3 = result3.T() * (-1);
    EXPECT_EQ(e1.z, result3);
}

TEST_F(MatrixS_TestCase, cat)
{
    MatrixS<2, 2> x = {1, 2, 3, 4};
    MatrixS<2, 2> y = {5, 6, 7, 8};
    MatrixS<2, 4> expect1{1, 2, 3, 4, 5, 6, 7, 8};
    MatrixS<4, 2> expect2{1, 2, 5, 6, 3, 4, 7, 8};
    EXPECT_EQ(expect1, concat<Ori::Col>(x, y));
    EXPECT_EQ(expect2, concat<Ori::Row>(x, y));

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

#ifdef PPX_USE_AVX

template <size_t M, size_t N, size_t L>
MatrixS<M, L> matmul(const MatrixS<M, N> &self, const MatrixS<N, L> &other)
{
    MatrixS<M, L> result;
    for (size_t k = 0; k < N; k++)
    {
        for (size_t j = 0; j < L; j++)
        {
            for (size_t i = 0; i < M; i++)
            {
                result(i, j) += self(i, k) * other(k, j);
            }
        }
    }
    return result;
}

TEST_F(MatrixS_TestCase, avx_full_mul)
{
    MatrixS<4, 4> A, B;
    random(A, -1e3, 1e3);
    random(B, -1e3, 1e3);
    EXPECT_EQ(A * B, matmul(A, B));

    MatrixS<4, 3> C;
    MatrixS<3, 4> D;
    random(C, -1e3, 1e3);
    random(D, -1e3, 1e3);
    EXPECT_EQ(C * D, matmul(C, D));

    MatrixS<4, 2> E;
    MatrixS<2, 4> F;
    random(E, -1e3, 1e3);
    random(F, -1e3, 1e3);
    EXPECT_EQ(E * F, matmul(E, F));

    MatrixS<4, 1> G;
    MatrixS<1, 4> H;
    random(G, -1e3, 1e3);
    random(H, -1e3, 1e3);
    EXPECT_EQ(G * H, matmul(G, H));

    MatrixS<16, 7> I;
    MatrixS<7, 16> J;
    random(I, -1e3, 1e3);
    random(J, -1e3, 1e3);
    EXPECT_EQ(I * J, matmul(I, J));

    MatrixS<80, 33> K;
    MatrixS<33, 20> L;
    random(K, -1e3, 1e3);
    random(L, -1e3, 1e3);
    EXPECT_EQ(K * L, matmul(K, L));

    MatrixS<100, 104> M;
    MatrixS<104, 100> N;
    random(M, -1e3, 1e3);
    random(N, -1e3, 1e3);
    EXPECT_EQ(M * N, matmul(M, N));
}

TEST_F(MatrixS_TestCase, avx_mask_mul)
{
    MatrixS<3, 3> A, B;
    random(A, -1e3, 1e3);
    random(B, -1e3, 1e3);
    EXPECT_EQ(A * B, matmul(A, B));

    MatrixS<3, 4> C;
    MatrixS<4, 3> D;
    random(C, -1e3, 1e3);
    random(D, -1e3, 1e3);
    EXPECT_EQ(C * D, matmul(C, D));

    MatrixS<3, 2> E;
    MatrixS<2, 3> F;
    random(E, -1e3, 1e3);
    random(F, -1e3, 1e3);
    EXPECT_EQ(E * F, matmul(E, F));

    MatrixS<3, 1> G;
    MatrixS<1, 3> H;
    random(G, -1e3, 1e3);
    random(H, -1e3, 1e3);
    EXPECT_EQ(G * H, matmul(G, H));

    MatrixS<3, 7> I;
    MatrixS<7, 3> J;
    random(I, -1e3, 1e3);
    random(J, -1e3, 1e3);
    EXPECT_EQ(I * J, matmul(I, J));

    MatrixS<3, 33> K;
    MatrixS<33, 3> L;
    random(K, -1e3, 1e3);
    random(L, -1e3, 1e3);
    EXPECT_EQ(K * L, matmul(K, L));

    MatrixS<3, 104> M;
    MatrixS<104, 3> N;
    random(M, -1e3, 1e3);
    random(N, -1e3, 1e3);
    EXPECT_EQ(M * N, matmul(M, N));
}

#endif

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}