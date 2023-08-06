from collections import deque
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pytransform3d.plot_utils import Frame
from pytransform3d.plot_utils import make_3d_axis
from pytransform3d import rotations as pr
from pytransform3d import transformations as pytr


# import modern_robotics as mr


class VisualPose:
    """Store Transform and trajectory for plot.

    func: function which will be called to update transform.
    """

    def __init__(self, axis, func, traj_len: int = 50):
        self.__func = func
        self.__trans = self.__func(0)
        self.__frame = Frame(self.__trans)
        self.__frame.add_frame(axis)

        self.__px = deque([self.__trans[0, 3]], maxlen=traj_len)
        self.__py = deque([self.__trans[1, 3]], maxlen=traj_len)
        self.__pz = deque([self.__trans[2, 3]], maxlen=traj_len)

        self.__line = axis.plot(self.__px, self.__py, self.__pz, c="r", alpha=0.2)[0]

    @property
    def trans(self):
        return self.__trans

    def __iter__(self):
        return (
            (self.__px[idx], self.__py[idx], self.__pz[idx])
            for idx in range(len(self.__px))
        )

    def __str__(self) -> str:
        return str(tuple(self))

    def update(self, timepoint: int, total_time: int):
        """update transform by single parameter."""
        self.__trans = self.__func(timepoint / total_time)
        self.__frame.set_data(self.__trans)
        self.__px.append(self.__trans[0, 3])
        self.__py.append(self.__trans[1, 3])
        self.__pz.append(self.__trans[2, 3])
        self.__line.set_data(self.__px, self.__py)
        self.__line.set_3d_properties(self.__pz)
        return self.__trans

    def clear(self):
        """clear trajectory."""
        self.__trans = self.__func(0)
        self.__frame.set_data(self.__trans)
        for pos in (self.__px, self.__py, self.__pz):
            pos.clear()
        return self.__trans


class VisualLink:
    """plot linkage between frame"""

    def __init__(self, axis, obj1: VisualPose, obj2: VisualPose):
        self.__line = axis.plot(
            (obj1.trans[0, 3], obj2.trans[0, 3]),
            (obj1.trans[1, 3], obj2.trans[1, 3]),
            (obj1.trans[2, 3], obj2.trans[2, 3]),
            c="k",
        )[0]
        self.__objs = obj1
        self.__obje = obj2

    def update(self):
        """update linkage with transform."""
        self.__line.set_data(
            (self.__objs.trans[0, 3], self.__obje.trans[0, 3]),
            (self.__objs.trans[1, 3], self.__obje.trans[1, 3]),
        )
        self.__line.set_3d_properties(
            (self.__objs.trans[2, 3], self.__obje.trans[2, 3])
        )


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
    transform_obj = []
    transform_obj.append(VisualPose(ax, cal_cur_trans1))
    transform_obj.append(VisualPose(ax, cal_cur_trans2))

    linkage_obj = []
    linkage_obj.append(VisualLink(ax, transform_obj[0], transform_obj[1]))

    def update_obj(cur_f):
        """Collect all animated objs to update."""
        for obj in transform_obj:
            obj.update(cur_f, NFRAME)
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
