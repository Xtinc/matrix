#include "gtest/gtest.h"
#include "LieGroup.hpp"

using namespace details;

class LieGroup : public testing::Test
{
};

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST_F(LieGroup, VecToSO3)
{
    Matrix<3, 1> vec{1, 2, 3};
    SO3 result{0, 3, -2, -3, 0, 1, 2, -1, 0};
    EXPECT_EQ(result, hat(vec));
}

TEST_F(LieGroup, MatrixLog6)
{
    SE3 Tinput{1, 0, 0, 0,
               0, 0, 1, 0,
               0, -1, 0, 0,
               0, 0, 3, 1};
    Matrix<4, 4> result{0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 1.57079633, 0.0,
                        0.0, -1.57079633, 0.0, 0.0,
                        0.0, 2.35619449, 2.35619449, 0.0};
    EXPECT_EQ(hat(Tinput.log()), result);
}

TEST_F(LieGroup, SE3Adt)
{
    SE3 T{1, 0, 0, 0,
          0, 0, 1, 0,
          0, -1, 0, 0,
          0, 0, 3, 1};
    Matrix<6, 6> result{1, 0, 0, 0, 3, 0,
                        0, 0, 1, 0, 0, 0,
                        0, -1, 0, 3, 0, 0,
                        0, 0, 0, 1, 0, 0,
                        0, 0, 0, 0, 0, 1,
                        0, 0, 0, 0, -1, 0};
    ASSERT_TRUE(result == T.Adt());
}

TEST_F(LieGroup, MatrixExp6)
{
    Matrix<4, 4> se3mat = {0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 1.5708, 0.0,
                           0.0, -1.5708, 0.0, 0.0,
                           0.0, 2.3562, 2.3562, 0.0};
    SE3 result{{1, 0, 0, 0, 0, 1, 0, -1, 0, 0}, {0, 0, 3}};
    auto cal = vee(se3mat).exp();
    for (size_t i = 0; i < 16; i++)
    {
        EXPECT_NEAR(cal[i], result[i], 1.0e-5);
    }
}