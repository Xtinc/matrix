import datetime
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pytransform3d.plot_utils import make_3d_axis
import modern_robotics as mr
from tools.misc import VisualPose
from tools.misc import VisualLink


class RobotTimeSeries:
    """struct to record robot operation"""

    def __init__(self, time, inst_q, real_q):
        self.time = time
        self.instq = inst_q
        self.realq = real_q

    def __iter__(self):
        return (c for c in (self.time, self.instq, self.realq))

    def __str__(self):
        return f"{self.time.isoformat()} inst:{self.instq} real:{self.realq}"


def convert_rts(b_str):
    """convert RobotTimeSeries object form string"""
    content = b_str.decode("UTF-8")
    time_point = datetime.datetime.fromisoformat(content[:26])
    idx1 = content.find("inst:[")
    idx2 = content.find("] real:[")
    inst_q = [float(qi) for qi in content[idx1 + 6 : idx2].split(",") if qi is not None]
    real_q = [float(qi) for qi in content[idx2 + 8 : -1].split(",") if qi is not None]
    return RobotTimeSeries(time_point, inst_q, real_q)


def register_rts():
    """init sqlite3 for RobotTimeSeries"""
    # Register the adapter and converter
    sqlite3.register_adapter(RobotTimeSeries, RobotTimeSeries.__str__)
    sqlite3.register_converter("RobotTimeSeries", convert_rts)


def init_database(filename: str):
    """Create a default SQLite3 db file."""
    connection = sqlite3.connect(filename, detect_types=sqlite3.PARSE_DECLTYPES)
    cursor = connection.execute("CREATE TABLE test(p RobotTimeSeries)")
    return connection, cursor


def close_database(connection):
    connection.commit()
    connection.close()


def insert_many(connection, rts_list):
    """insert a list of RobotTimeSeries objs"""
    # Successful, con.commit() is called automatically afterwards
    with connection:
        connection.executemany("INSERT INTO test(p) VALUES(?)", rts_list)


def insert_one(connection, rts):
    """insert RobotTimeSeries objs"""
    connection.execute("INSERT INTO test(p) VALUES(?)", (rts,))


def read_database(connection):
    """read list of RobotTimeSeries objs from database"""
    return [
        row for row in connection.execute('SELECT p AS "p [RobotTimeSeries]" FROM test')
    ]


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


def ur6_joint_gen(timestep):
    """generator of ur6 joint list"""
    return np.ones(6) * np.pi * 2 * timestep


if __name__ == "__main__":
    plt.style.use("tools.gnuplot")
    fig = plt.figure()
    ax = make_3d_axis(ax_s=800)

    NFRAME = 100
    register_rts()
    con = sqlite3.connect("test.sqlite", detect_types=sqlite3.PARSE_DECLTYPES)
    con.execute("CREATE TABLE test(p RobotTimeSeries)")

    tlist = fk_space_all(np.array([0, 0, 0, 0, 0, 0]))
    transform_obj = [VisualPose(ax, t, traj_len=0, scale=30) for t in tlist]
    linkage_obj = []
    for i in range(len(tlist) - 1):
        linkage_obj.append(VisualLink(ax, transform_obj[i], transform_obj[i + 1]))

    def update_obj(cur_f):
        """Collect all animated objs to update."""
        theta_list = ur6_joint_gen(cur_f / NFRAME)
        insert_one(
            con, RobotTimeSeries(datetime.datetime.now(), theta_list, theta_list)
        )
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
    close_database(con)
