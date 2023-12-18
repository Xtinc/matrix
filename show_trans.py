"""Simulated UR6 Robotics"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pytransform3d.plot_utils import make_3d_axis
import modern_robotics as mr
from tools.misc import VisualPose
from tools.misc import VisualLink


if __name__ == "__main__":
    plt.style.use("tools.gnuplot")
    fig = plt.figure()
    ax = make_3d_axis(ax_s=2)

    NFRAME = 200
    trans_strat = mr.MatrixExp6(mr.VecTose3([0, 0, 0, 0, 0, 1]))
    trans_end = mr.MatrixExp6(mr.VecTose3([np.pi / 2, np.pi / 2, np.pi / 2, 0, 0, 0]))

    transform_start = VisualPose(ax, trans_strat, traj_len=0, scale=0.5)
    transform_end = VisualPose(ax, trans_end, traj_len=0, scale=0.5)
    obj = VisualPose(ax, trans_strat, traj_len=50, scale=0.5)

    def update_obj(cur_f):
        """Collect all animated objs to update."""
        screw = mr.MatrixLog6(mr.TransInv(trans_strat) @ trans_end)
        vec_screw = mr.se3ToVec(screw)
        vec_screw[1] = 0.0
        screw = mr.VecTose3(vec_screw)
        tend = mr.MatrixExp6(cur_f / NFRAME * screw)
        obj.update(trans_strat @ tend)

    def init_obj():
        """Collect all animated objs to init."""
        obj.clear()
        obj.update(trans_strat)

    anim = animation.FuncAnimation(
        fig,
        update_obj,
        NFRAME,
        init_obj,
        interval=25,
    )

    plt.show()
