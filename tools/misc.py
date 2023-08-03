"""A cluster of miscellaneous things.

Contains:
    Packed FFT function.
"""
import numpy as np
from scipy.fft import fft as scipyfft
from scipy.signal import welch as scipywelch


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


def convertq(content):
    """Convert Q from a string, size should be same as q."""
    list_q = [float(ele) for ele in content.split()]
    return list_q
