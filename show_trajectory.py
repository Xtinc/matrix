import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pytransform3d.plot_utils import make_3d_axis
from pytransform3d import rotations as pr
from pytransform3d import transformations as pytr
import modern_robotics as mr
from tools.misc import VisualPose
from tools.misc import VisualLink

R1 = pr.active_matrix_from_extrinsic_euler_zyz([np.pi / 3, np.pi / 6, np.pi / 2])
R2 = pr.active_matrix_from_extrinsic_euler_zyz([-np.pi / 3, np.pi / 6, -np.pi / 2])

T1 = pytr.transform_from(R=R1, p=[1, 1, 1])
T2 = pytr.transform_from(R=R2, p=[-1, -1, -1])
T3 = pytr.transform_from(R=R2, p=[-1.5, -1.5, -1.5])
dis = pytr.exponential_coordinates_from_transform_log(
    pytr.transform_log_from_transform(np.matmul(T1, np.linalg.inv(T2)))
)


def cal_cur_trans1(timestep):
    """Calculate Transform by matrix exp"""
    distance = pytr.transform_from_exponential_coordinates(timestep * dis)
    return np.matmul(distance, T2)


def cal_cur_trans2(timestep):
    """Calculate Transform by matrix exp"""
    distance = pytr.transform_from_exponential_coordinates(timestep * dis)
    return np.matmul(distance, T3)


if __name__ == "__main__":
    plt.style.use("tools.gnuplot")
    fig = plt.figure(figsize=(5, 5))
    ax = make_3d_axis(ax_s=2)

    pytr.plot_transform(ax, A2B=T2, name="Start")
    pytr.plot_transform(ax, A2B=T1, name="End")

    NFRAME = 200

    function_list = []
    function_list.append(cal_cur_trans1)
    function_list.append(cal_cur_trans2)

    transform_obj = []
    transform_obj.append(VisualPose(ax, cal_cur_trans1(0)))
    transform_obj.append(VisualPose(ax, cal_cur_trans2(0)))

    linkage_obj = []
    linkage_obj.append(VisualLink(ax, transform_obj[0], transform_obj[1]))

    def update_obj(cur_f):
        """Collect all animated objs to update."""
        for func, obj in zip(function_list, transform_obj):
            obj.update(func(cur_f / NFRAME))
        for link in linkage_obj:
            link.update()

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
        interval=int(10000 / NFRAME),
    )

    plt.show()
