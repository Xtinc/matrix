#ifndef VVERY_SIMPLE_ROBOTICS_HEADER
#define VVERY_SIMPLE_ROBOTICS_HEADER

#include "liegroup.hpp"

namespace ppx
{
    using range = std::pair<double, double>;
    struct joint
    {
        std::string name{"Anon"};
        range range{-gl_rep_max, gl_rep_max};
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
        using idx_available_t = typename std::enable_if<gl_less_than(L, N), RT>::type;
        std::array<joint, N> JList;

    public:
        using Q = Matrix<N, 1>;
        template <size_t L>
        idx_available_t<L> setJoint(const joint &j)
        {
            JList[L] = j;
        }
        SE3 forwardSpace(const std::string &jointName, const std::vector<double> &jointAngle) const
        {
            // check name in joints list & calculate poe.
            int real_size = (int)gl_get_less_dynamic(N, jointAngle.size());
            int effector_idx = -1;
            for (int i = 0; i < real_size; i++)
            {
                if (JList[i].name == jointName)
                {
                    effector_idx = i;
                    break;
                }
            }
            if (effector_idx == -1)
            {
                return {};
            }
            SE3 effector_pose = JList[effector_idx].pose;
            for (int i = effector_idx; i > -1; i--)
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
        // Matrix<6, N> jacobiSpace(const std::vector<double> &jointAngle)
        // {
        //     int real_size = (int)gl_get_less_dynamic(N, jointAngle.size());
        //     std::array<double, N> JAngle{};
        //     std::copy_n(jointAngle.begin(), real_size, JAngle.begin());
        //     Matrix<6, N> Js;
        //     SE3 T;
        //     for (int i = 1; i < real_size; i++)
        //     {
        //         se3 screw = JList[i - 1].screw * JAngle[i - 1];
        //         T = T * screw.exp();
        //         Js({-1, -1}, i) = T.Adt() * JList[i].screw;
        //     }
        //     return Js;
        // }
        Matrix<6, N> jacobiSpace(const Q &jointAngle)
        {
            Matrix<6, N> Js;
            SE3 T;
            for (int i = 1; i < N; i++)
            {
                se3 screw = JList[i - 1].screw * jointAngle[i - 1];
                T = T * screw.exp();
                Js({-1, -1}, i) = T.Adt() * JList[i].screw;
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
                auto err_w = norm2(Vs.w());
                auto err_v = norm2(Vs.v());
                // printf("iter=%d, w_error=%f, v_error=%f\n", iter, err_w, err_v);
                converage = err_w < gl_rep_eps && err_v < gl_rep_eps;
                init += pinv(jacobiSpace(init)) * Vs;
                ++iter;
            }
            return init;
        }
    };
}

#endif