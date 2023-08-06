"""A cluster of miscellaneous things.

Contains:
    Packed FFT function.
"""
from collections import deque
import numpy as np
from scipy.fft import fft as scipyfft
from scipy.signal import welch as scipywelch
from pytransform3d.plot_utils import Frame


def fft(x_s, f_s):
    """Calculate FFT of sequence xs, fs is sampling frequence."""
    idx = int(np.log2(len(x_s)))
    seq_x = x_s[1 : int(np.exp2(idx))]
    data_len = len(seq_x)
    half_len = int(data_len / 2)
    fft_x = np.array(scipyfft(seq_x, norm="forward"))
    amp_x = abs(fft_x) * 2
    lable_x = np.linspace(0, half_len - 1, half_len)
    amp_half = amp_x[0:half_len]
    amp_half[0] = 0
    fre = lable_x / data_len * f_s
    pha = np.unwrap(np.angle(fft_x))
    return fre, amp_half, pha


def psd(x_s, f_s, scale_type="density"):
    """Calculate PSD by welch Periodogram, fs is sampling frequence."""
    idx = int(np.log2(len(x_s)))
    if scale_type == "density":
        return scipywelch(x_s, f_s, nperseg=int(np.exp2(idx)), scaling="density")
    return scipywelch(
        x_s, f_s, "flattop", nperseg=int(np.exp2(idx)), scaling="spectrum"
    )


class VisualPose:
    """Store Transform and trajectory for plot.

    func: function which will be called to update transform.
    """

    def __init__(self, axis3d, trans, traj_len: int = 50, scale=1.0):
        self.__trans = trans
        self.__frame = Frame(self.__trans, s=scale)
        self.__frame.add_frame(axis3d)

        self.__px = deque([self.__trans[0, 3]], maxlen=traj_len)
        self.__py = deque([self.__trans[1, 3]], maxlen=traj_len)
        self.__pz = deque([self.__trans[2, 3]], maxlen=traj_len)

        self.__line = axis3d.plot(self.__px, self.__py, self.__pz, c="r", alpha=0.2)[0]

    @property
    def trans(self):
        """return transform now."""
        return self.__trans

    def __iter__(self):
        return (
            (self.__px[idx], self.__py[idx], self.__pz[idx])
            for idx in range(len(self.__px))
        )

    def __str__(self) -> str:
        return str(tuple(self))

    def update(self, trans):
        """update transform by single parameter."""
        self.__trans = trans
        self.__frame.set_data(self.__trans)
        self.__px.append(self.__trans[0, 3])
        self.__py.append(self.__trans[1, 3])
        self.__pz.append(self.__trans[2, 3])
        self.__line.set_data(self.__px, self.__py)
        self.__line.set_3d_properties(self.__pz)

    def clear(self):
        """clear trajectory."""
        self.__frame.set_data(self.__trans)
        for pos in (self.__px, self.__py, self.__pz):
            pos.clear()


class VisualLink:
    """plot linkage between frame"""

    def __init__(self, axis3d, obj1: VisualPose, obj2: VisualPose):
        self.__line = axis3d.plot(
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
