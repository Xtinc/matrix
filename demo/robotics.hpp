#ifndef VVERY_SIMPLE_ROBOTICS_HEADER
#define VVERY_SIMPLE_ROBOTICS_HEADER

#include "liegroup.hpp"

using namespace ppx;

struct joint
{
    std::string name{"Anon"};
    std::pair<double, double> range{-MAX_SP, MAX_SP};
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

#endif