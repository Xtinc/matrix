#include "robotics.h"

using namespace ppx;

namespace robotics
{
    SE3 kinematics::forwardSpace(const std::string &jointName, const std::vector<double> &jointAngle) const
    {
        // check name in joints list & calculate poe.
        auto real_size = gl_get_less_dynamic(jointLists.size(), jointAngle.size());
        auto effector_idx = -1;
        for (size_t i = 0; i < real_size; i++)
        {
            if (jointLists[i].jointName == jointName)
            {
                effector_idx = static_cast<int>(i);
                break;
            }
        }
        if (effector_idx == -1)
        {
            return {};
        }
        SE3 effector_pose = jointLists[effector_idx].initPose;
        for (int i = effector_idx; i > -1; i--)
        {
            se3 screw = jointLists[i].jointScrew * jointAngle[i];
            effector_pose = screw.exp() * effector_pose;
        }
        return effector_pose;
    }

}