#ifndef VVERY_SIMPLE_ROBOTICS_HEADER
#define VVERY_SIMPLE_ROBOTICS_HEADER

#include "liegroup.hpp"

namespace robotics
{
    using range = std::pair<double, double>;
    struct joint
    {
        std::string jointName{"Anon"};
        range jointRange{-ppx::gl_rep_max, ppx::gl_rep_max};
        ppx::se3 jointScrew{};
        ppx::SE3 initPose{};
        joint() = default;
        joint(const std::string &name, const ppx::se3 &screw, const ppx::SE3 &pose)
            : jointName(name), jointScrew(screw), initPose(pose){};
    };

    class kinematics
    {
    private:
        std::vector<joint> jointLists;

    public:
        void pushJoint(const joint &j)
        {
            jointLists.push_back(j);
        }
        // ppx::SE3 forwardBody(const std::string &jointName, const std::vector<double> &jointAngle) const;
        ppx::SE3 forwardSpace(const std::string &jointName, const std::vector<double> &jointAngle) const;
    };
}

#endif