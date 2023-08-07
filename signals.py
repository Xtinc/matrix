import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import tools.misc as misc

plt.style.use("tools.gnuplot")
N = 2048 * 5  # 采样点的个数
x = np.arange(0, 2 * np.pi, 2 * np.pi / N)
fs = N / (2 * np.pi)
# 产生频率为120、500、10hz的信号进行模拟
y = (
    7 * np.sin(120 * 2 * np.pi * x)
    + 5 * np.sin(500 * 2 * np.pi * x)
    + 9 * np.sin(10 * 2 * np.pi * x)
    + np.random.normal(scale=np.sqrt(1.89), size=len(x))
)
fre, amp, _ = misc.fft(y, fs)
plt.plot(fre, amp)
# w = np.arange(0, N, 1)  # 频域轴


# b1 = signal.firwin(
#     51, [0.42, 0.8], window="hamming", pass_zero="bandstop"
# )  # 哈明窗，截至频率100Hz
# b2 = [
#     0.0010176,
#     0.000927711,
#     -0.00238796,
#     0.000434621,
#     -0.00010929,
#     0.00253823,
#     0.00180104,
#     -0.00843451,
#     0.00289083,
#     0.00124462,
#     0.00682364,
#     0.00227208,
#     -0.0243103,
#     0.0124667,
#     0.00617722,
#     0.0127504,
#     -0.00127084,
#     -0.0580584,
#     0.0431386,
#     0.0189968,
#     0.0179269,
#     -0.0191777,
#     -0.172457,
#     0.224718,
#     0.120414,
#     0.619334,
#     0.120414,
#     0.224718,
#     -0.172457,
#     -0.0191777,
#     0.0179269,
#     0.0189968,
#     0.0431386,
#     -0.0580584,
#     -0.00127084,
#     0.0127504,
#     0.00617722,
#     0.0124667,
#     -0.0243103,
#     0.00227208,
#     0.00682364,
#     0.00124462,
#     0.00289083,
#     -0.00843451,
#     0.00180104,
#     0.00253823,
#     -0.00010929,
#     0.000434621,
#     -0.00238796,
#     0.000927711,
#     0.0010176,
# ]
# w1, h = signal.freqz(b1)  # 求频响
# w2, h2 = signal.freqz(b2)

# plt.figure(1)
# plt.title("Frequence Response")
# plt.plot(w1 / 2 / np.pi * N, 20 * np.log10(np.abs(h) + 0.01))
# plt.xlabel("$f$")
# plt.ylabel("$H(f)$")
# plt.plot(w2 / 2 / np.pi * N, 20 * np.log10(np.abs(h2) + 0.01))

# b2 = signal.firwin(24, 2 * 100 / N, window="hann")  # 汉宁窗，截至频率100Hz
# w1, h = signal.freqz(b2)  # 求频响
# plt.figure(2)
# plt.title("freqz")
# plt.plot(w1 / 2 / np.pi * N, 20 * np.log10(np.abs(h) + 0.01))

# b3 = signal.firwin(24, 2 * 100 / N, window="blackman")  # 布莱克曼窗，截至频率100Hz
# w1, h = signal.freqz(b3)  # 求频响
# plt.figure(3)
# plt.title("freqz")
# plt.plot(w1 / 2 / np.pi * N, 20 * np.log10(np.abs(h) + 0.01))

# b4 = signal.firwin(N, 2 * 100 / N, window="boxcar")  # 矩形窗，截至频率100Hz
# w1, h = signal.freqz(b4)  # 求频响
# plt.figure(4)
# plt.title("freqz")
# plt.plot(w1 / 2 / np.pi * N, 20 * np.log10(np.abs(h) + 0.01))
plt.show()
