#include "gtest/gtest.h"
#include "liegroup.hpp"

using namespace ppx;

class MatrixS_TestCase : public ::testing::Test
{
public:
    MatrixS_TestCase(){};
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
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}