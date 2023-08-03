"""A cluster of miscellaneous things.

Contains:
    Packed FFT function.
"""
import numpy as np
from scipy.fft import fft as scipyfft


def fft(x_s, f_s):
    """Calculate FFT of sequence xs , fs is sampling frequence."""
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
    return amp_half, fre, pha


def signal_gen(sample_num, sample_freq, signal_type):
    """generator of a signal sequence. type should be uniform/guassian"""
    signal_sequence = []
    if signal_type == "linear":
        signal_sequence = np.arange(sample_num) / sample_freq
    elif signal_type == "guassian":
        signal_sequence = np.random.normal(size=sample_num)
    return signal_sequence


def convertq(content):
    """Convert Q from a string, size should be same as q."""
    list_q = [float(ele) for ele in content.split()]
    return list_q
