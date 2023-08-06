import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pytransform3d.plot_utils import make_3d_axis
import modern_robotics as mr
from tools.misc import VisualPose
from tools.misc import VisualLink


M01 = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 89.159], [0, 0, 0, 1]]
M12 = [[0, 0, 1, 0], [0, 1, 0, 135.85], [-1, 0, 0, 0], [0, 0, 0, 1]]
M23 = [[1, 0, 0, 0], [0, 1, 0, -119.7], [0, 0, 1, 425], [0, 0, 0, 1]]
M34 = [[0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 392.25], [0, 0, 0, 1]]
M45 = [[1, 0, 0, 0], [0, 1, 0, 93], [0, 0, 1, 0], [0, 0, 0, 1]]
M56 = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 94.65], [0, 0, 0, 1]]
M67 = [[1, 0, 0, 0], [0, 0, 1, 82.3], [0, -1, 0, 0], [0, 0, 0, 1]]


Mlist_R = [M01, M12, M23, M34, M45, M56, M67]
Mlist_G = [np.eye(4, 4)]
ME = Mlist_G[0]
for ele in Mlist_R:
    ME = np.matmul(ME, ele)
    Mlist_G.append(ME)
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


def fk_space_all(thetalist):
    """Calculate UR6 FK all."""
    result = []
    for idx, endeffector in enumerate(Mlist_G):
        if idx == 0:
            result.append(Mlist_G[0])
        else:
            result.append(mr.FKinSpace(endeffector, Slist[:, :idx], thetalist[:idx]))
    return result


if __name__ == "__main__":
    plt.style.use("tools.gnuplot")
    fig = plt.figure()
    ax = make_3d_axis(ax_s=800)

    NFRAME = 100

    tlist = fk_space_all(np.array([0, 0, 0, 0, 0, 0]))
    transform_obj = [VisualPose(ax, t, traj_len=0, scale=30) for t in tlist]
    linkage_obj = []
    for i in range(len(tlist) - 1):
        linkage_obj.append(VisualLink(ax, transform_obj[i], transform_obj[i + 1]))

    def update_obj(cur_f):
        """Collect all animated objs to update."""
        theta_list = np.ones(6) * np.pi * 2 * cur_f / NFRAME
        tlist_cur = fk_space_all(theta_list)
        for tend, obj in zip(tlist_cur, transform_obj):
            obj.update(tend)
        for link in linkage_obj:
            link.update()

    def init_obj():
        """Collect all animated objs to init."""
        for obj in transform_obj:
            obj.clear()
        for link in linkage_obj:
            link.update()

    anim = animation.FuncAnimation(
        fig, update_obj, NFRAME, init_obj, interval=int(10000 / NFRAME)
    )

    plt.show()
