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
    Q inverseSpace(const SE3 &pose, Q init)
    {
        auto Tsb = forwardSpace(init);
        auto Js = jacobiSpace(init);
        auto err = (pose * Tsb.I()).log();

        auto dtd = 4 * inner_product(err, err);
        auto iter = 0u;

        while (++iter < 50 && norm2(err) > EPS_SP && sqrt(dtd) > 100 * EPS_SP)
        {
            auto err_w = norm2(err._1());
            auto err_v = norm2(err._2());
            printf("iter=%d, w_error=%15.8f, v_error=%15.8f, err=%15.8f\n", iter, err_w, err_v, norm2(err));

            auto g = Js.T() * err;
            auto alpha = inner_product(err, Js * g) / SQR(norm2(Js * g));
            auto result = linsolve<Factorization::LU>(Js.T() * Js, g);

            Q dq{};
            if (result.s == StatusCode::CONVERGED)
            {
                Q pU = alpha * g;
                const auto &pB = result.x;
                auto npB = norm2(pB);
                auto npU = norm2(pU);
                if (npB < dtd)
                {
                    std::cout << "use GN ";
                    dq = result.x;
                }
                else if (npU > dtd)
                {
                    std::cout << "use GD1 ";
                    dq = sqrt(dtd) / norm2(g) * g;
                }
                else
                {
                    std::cout << "use MX ";
                    Q pBU = pB - pU;
                    auto tau = sqrt(SQR(inner_product(pU, pBU)) - inner_product(pBU, pBU) * (npU - dtd));
                    tau = (tau - inner_product(pU, pBU)) / inner_product(pBU, pBU);
                    dq = pU + tau * pBU;
                }
            }
            else
            {
                std::cout << "use GD2 ";
                dq = sqrt(dtd) / norm2(g) * g;
                // should use GD in trust region ?
            }

            auto pseudo_init = init + dq;
            auto pseudo_Tsb = forwardSpace(pseudo_init);
            auto pseudo_Js = jacobiSpace(pseudo_init);
            auto pseudo_err = (pose * pseudo_Tsb.I()).log();
            auto rho = inner_product(err, err) - inner_product(pseudo_err, pseudo_err);
            rho /= inner_product(dq, Js.T() * err);
            std::cout << "rho: " << rho << " dtd: " << dtd << "\n";
            if (rho > 1e-5)
            {
                init = pseudo_init;
                Tsb = pseudo_Tsb;
                Js = pseudo_Js;
                err = pseudo_err;
            }
            if (rho > 0.75)
            {
                dtd = std::max(dtd, 9 * norm2(dq));
            }
            else if (rho < 0.25)
            {
                dtd /= 9;
            }
        }
        return init;
    }
};

#endif