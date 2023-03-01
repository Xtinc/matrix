#include "gtest/gtest.h"
#include "liegroup.hpp"

using namespace ppx;

struct joint
{
    std::string name{"Anon"};
    se3 screw;
    SE3 pose;
    joint() = default;
    joint(const std::string &name_, const se3 &screw_, const SE3 &pose_)
        : name(name_), screw(screw_), pose(pose_){};
};

template <size_t N>
class kinematics
{
private:
    template <size_t L, typename RT = void>
    using idx_available_t = typename std::enable_if<(L < N), RT>::type;
    std::array<joint, N> JList;

public:
    using Q = MatrixS<N, 1>;
    template <size_t L>
    idx_available_t<L, joint &> getJoint()
    {
        return JList[L];
    }
    template <size_t L>
    idx_available_t<L, const joint &> getJoint() const
    {
        return JList[L];
    }
    template <size_t L>
    idx_available_t<L> setJoint(const joint &j)
    {
        JList[L] = j;
    }
    SE3 forwardSpace(const std::string &jointName, const Q &jointAngle) const
    {
        auto effector_idx = std::find_if(JList.begin(), JList.end(), [&jointName](const joint &elem)
                                         { return jointName == elem.name; });
        if (effector_idx == JList.end())
        {
            return {};
        }
        SE3 effector_pose = effector_idx->pose;
        for (int i = (int)std::distance(JList.begin(), effector_idx); i > -1; i--)
        {
            se3 screw = JList[i].screw * jointAngle[i];
            effector_pose = screw.exp() * effector_pose;
        }
        return effector_pose;
    }
    SE3 forwardSpace(const Q &jointAngle) const
    {
        SE3 effector_pose = JList.back().pose;
        for (int i = N - 1; i > -1; i--)
        {
            se3 screw = JList[i].screw * jointAngle[i];
            effector_pose = screw.exp() * effector_pose;
        }
        return effector_pose;
    }
    MatrixS<6, N> jacobiSpace(const Q &jointAngle)
    {
        MatrixS<6, N> Js;
        SE3 T;
        for (int i = 1; i < (int)N; i++)
        {
            se3 screw = JList[i - 1].screw * jointAngle[i - 1];
            T = T * screw.exp();
            Js.template sub<6, 1>(0, i) = T.Adt() * JList[i].screw;
        }
        return Js;
    }
    template <size_t L>
    idx_available_t<L, MatrixS<6, L>> jacobiSpace(const std::array<std::string, L> &namelist, const Q &jointAngle)
    {
        MatrixS<6, L> Js;
        auto JsAll = jacobiSpace(jointAngle);
        for (size_t i = 0; i < L; i++)
        {
            auto iter = std::find_if(JList.begin(), JList.end(), [&namelist, i](const joint &elem)
                                     { return namelist[i] == elem.name; });
            if (iter != JList.end())
            {
                auto col_idx = std::distance(JList.begin(), iter);
                Js.template sub<6, 1>(0, i) = JsAll.template sub<6, 1>(0, col_idx);
            }
        }
        return Js;
    }
    Q inverseSpace(const SE3 &pose, Q init)
    {
        SE3 Tsb;
        se3 Vs;
        bool converage = false;
        auto iter = 0u;
        while (!converage && iter < 20)
        {
            Tsb = forwardSpace(init);
            Vs = Tsb.Adt() * (Tsb.I() * pose).log();
            auto err_w = norm2(Vs._1());
            auto err_v = norm2(Vs._2());
            printf("iter=%d, w_error=%f, v_error=%f\n", iter, err_w, err_v);
            converage = err_w < EPS_SP && err_v < EPS_SP;
            init += linsolve<Factorization::SVD>(jacobiSpace(init), Vs).x;
            ++iter;
        }
        return init;
        // auto fn = [this, pose](const Q &x)
        // {
        //     auto Tsb = forwardSpace(x);
        //     auto Vs = Tsb.Adt() * (Tsb.I() * pose).log();
        //     return 0.5 * inner_product(Vs, Vs);
        // };
        // auto dfn = [this, pose](const Q &x)
        // {
        //     auto Tsb = forwardSpace(x);
        //     se3 Vs = Tsb.Adt() * (Tsb.I() * pose).log();
        //     return Q(jacobiSpace(x) * Vs);
        // };
        // return fminunc<Optimization::GradientDescent>(fn, dfn, init).x;
    }
};

class LieGroup_TestCase : public ::testing::Test
{
public:
    LieGroup_TestCase() = default;
};

TEST_F(LieGroup_TestCase, so3)
{
    MatrixS<3, 1> vec{1, 2, 3};
    SO3 expected{0, 3, -2, -3, 0, 1, 2, -1, 0};
    SO3 result = hat(vec);
    EXPECT_EQ(result, expected);
}

TEST_F(LieGroup_TestCase, SE3)
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
}

TEST_F(LieGroup_TestCase, se3)
{
    MatrixS<4, 4> se3mat = {0.0, 0.0, 0.0, 0.0,
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

TEST_F(LieGroup_TestCase, kinematics)
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
    SE3 result{0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.095, 0.109, 0.988, 1.0};
    auto r = UR5.forwardSpace("R6", {0.0, -0.5 * PI, 0.0, 0.0, 0.5 * PI, 0.0});
    EXPECT_EQ(r, result);
    MatrixS<6, 6> result2{0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 1.0, 0.0, -0.089, 0.0, 0.0,
                          0.0, 1.0, 0.0, -0.514, 0.0, 0.0,
                          0.0, 1.0, 0.0, -0.906, 0.0, 0.0,
                          1.0, 0.0, 0.0, 0.0, 0.906, -0.109,
                          0.0, 0.0, 1.0, 0.109, -0.095, 0.0};
    auto r2 = UR5.jacobiSpace({0.0, -0.5 * PI, 0.0, 0.0, 0.5 * PI, 0.0});
    EXPECT_EQ(r2, result2);
    auto r3 = UR5.jacobiSpace(std::array<std::string, 3>{"R1", "R2", "R3"}, {0.0, -0.5 * PI, 0.0, 0.0, 0.5 * PI, 0.0});
    MatrixS<6, 3> result3 = r2.sub<6, 3>(0, 0);
    EXPECT_EQ(r3, result3);
    SE3 TargetPose{0.0, 1.0, 0.0, 0.0,
                   -1.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 1.0, 0.0,
                   0.095, 0.109, 0.988, 1.0};
    auto j = UR5.inverseSpace(TargetPose, {0.0, -1.5, 0.0, 0.0, 1.5, 0.0});
    EXPECT_EQ(UR5.forwardSpace("R6", j), TargetPose);
}