"""Simulated UR6 Robotics"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pytransform3d.plot_utils import make_3d_axis
import modern_robotics as mr
from tools.misc import VisualPose
from tools.misc import VisualLink

M01 = np.array(
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 89.159],
        [0.0, 0.0, 0.0, 1.0],
    ]
)
M02 = np.array(
    [
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 135.85],
        [-1.0, 0.0, 0.0, 89.159],
        [0.0, 0.0, 0.0, 1.0],
    ]
)
M03 = np.array(
    [
        [0.0, 0.0, 1.0, 425.0],
        [0.0, 1.0, 0.0, 16.15],
        [-1.0, 0.0, 0.0, 89.159],
        [0.0, 0.0, 0.0, 1.0],
    ]
)
M04 = np.array(
    [
        [-1.0, 0.0, 0.0, 817.25],
        [0.0, 1.0, 0.0, 16.15],
        [0.0, 0.0, -1.0, 89.159],
        [0.0, 0.0, 0.0, 1.0],
    ]
)
M05 = np.array(
    [
        [-1.0, 0.0, 0.0, 817.25],
        [0.0, 1.0, 0.0, 109.15],
        [0.0, 0.0, -1.0, 89.159],
        [0.0, 0.0, 0.0, 1.0],
    ]
)
M06 = np.array(
    [
        [-1.0, 0.0, 0.0, 817.25],
        [0.0, 1.0, 0.0, 109.15],
        [0.0, 0.0, -1.0, -5.491],
        [0.0, 0.0, 0.0, 1.0],
    ]
)


Mlist = [np.eye(4, 4), M01, M02, M03, M04, M05, M06]
Slist = np.array(
    [
        [0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 0, 1],
        [1, 0, 0, 0, -1, 0],
        [0, -89.159, -89.159, -89.159, -109.15, 0.5491],
        [0, 0, 0, 0, 817.25, 0],
        [0, 0, 425, 817.25, 0, 817.25],
    ]
)


def fk_space(thetalist):
    """Calculate UR6 FK all."""
    result = [np.eye(4, 4)]
    for idx, endeffector in enumerate(Mlist):
        if idx != 0:
            result.append(mr.FKinSpace(endeffector, Slist[:, :idx], thetalist[:idx]))
    return result


def ik_space(pose, thetalist0):
    """Inverse Kinematics"""
    return mr.IKinSpace(Slist, M06, pose, thetalist0, eomg=0.01, ev=0.001)


def ur6_joint_gen(timestep):
    """generator of ur6 joint list"""
    return np.ones(6) * np.pi * 2 * timestep


if __name__ == "__main__":
    plt.style.use("tools.gnuplot")
    fig = plt.figure()
    ax = make_3d_axis(ax_s=800)

    NFRAME = 200
    trans_list = []
    for idx in range(NFRAME):
        theta_list = ur6_joint_gen(idx / NFRAME)
        trans_list.append(fk_space(theta_list))

    transform_obj = [VisualPose(ax, t, traj_len=0, scale=0) for t in Mlist[:-1]]
    transform_obj.append(VisualPose(ax, Mlist[-1], traj_len=120, scale=100))
    linkage_obj = []
    for i in range(len(transform_obj) - 1):
        linkage_obj.append(VisualLink(ax, transform_obj[i], transform_obj[i + 1]))

    theta_list0 = np.ones(6)

    def update_obj(cur_f, theta0, _):
        """Collect all animated objs to update."""
        theta, sts = ik_space(cur_f, theta0)
        if not sts:
            print("Diverge")
        trans_list.append(fk_space(theta))
        for tend, obj in zip(trans_list[cur_f], transform_obj):
            obj.update(tend)
        for link in linkage_obj:
            link.update()
        theta0 = theta

    def init_obj():
        """Collect all animated objs to init."""
        for obj in transform_obj:
            obj.clear()
        for link in linkage_obj:
            link.update()

    anim = animation.FuncAnimation(
        fig,
        update_obj,
        NFRAME,
        init_obj,
        fargs=(theta_list0, NFRAME),
        interval=int(10000 / NFRAME),
    )

    plt.show()
