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
    joint(const std::string &name_, const se3 &screw_, const SE3 &pose_, double lo_, double hi_)
        : joint(name_, screw_, pose_)
    {
        range.first = lo_;
        range.second = hi_;
    };

    friend std::ostream &operator<<(std::ostream &os, const joint &self)
    {
        std::cout << "Name:\t" << self.name << "\nrange:\t" << self.range.first << " " << self.range.second
                  << "\nscrew:\t" << self.screw << "\npose:\t" << self.pose;
        return os;
    }
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
    idx_available_t<L, joint &> Joint()
    {
        return JList[L];
    }

    template <size_t L>
    idx_available_t<L, const joint &> Joint() const
    {
        return JList[L];
    }

    std::array<joint, N> &JointList()
    {
        return JList;
    }

    const std::array<joint, N> &JointList() const
    {
        return JList;
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
        Js.template sub<6, 1, false>(0, 0) = JList[0].screw;
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

    Q inverseSpace(const SE3 &pose, const Q &init)
    {
        auto fx = [this, &pose](const Q &q)
        {
            auto Tsb = this->forwardSpace(q);
            return (Tsb * pose.I()).log();
        };

        auto dfx = [this, &pose](const Q &q, const se3 &)
        {
            auto Jp = eye<6>();
            auto Tsb = this->forwardSpace(q) * pose.I();
            Jp.sub<3, 3, false>(3, 0) = hat(Tsb.Pos());
            return Jp.I() * this->jacobiSpace(q);
        };

        Q lower, upper;
        for (size_t i = 0; i < N; i++)
        {
            lower[i] = JList[i].range.first;
            upper[i] = JList[i].range.second;
        }

        lower.fill(-PI);
        upper.fill(PI);
        details::CoDo<Q::LEN, 6> codo(fx, dfx, lower, upper);

        auto result = codo(init);
        return result.s == StatusCode::CONVERGED ? result.x : init;
    }
};

#endif