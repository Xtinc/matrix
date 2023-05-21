# Chapter 4

## 频率成型滤波器的一般程序设计General program design of frequency-shaping filter



### 处理随机信号的两种途径

对于某随机过程，预测其随时间演化过程是常见的任务。通常会使用先验知识与模型然后基于概率推演，这也是一般随机过程分析常用的思路，如贝叶斯分析[^1]。但当随机过程中关注的变量间有较强的独立性与周期性时，也可以通过傅里叶变换在频率域中直接操作其频谱，根据需要进行具体的频率频率成型操作[^2]。

对于这类时间演化系统的估计，根据使用数据因果性，可具体分成三类任务：

- 滤波（filter）: 基于系统前n时刻的历史数据（含当前）估计下一时刻数据。

- 平滑（smooth）: 收集系统一段时间n内数据，根据所有数据，得到每个时刻的数据。

- 预测（predict）: 基于系统前n时刻的历史数据（含当前）估计n+t刻数。

三类任务本质上同种同源，只是根据时间因果性进行区分。因此很多算法：如Savizkg-Golay平滑算法被称作Savizkg-Golay滤波器；均值滤波器的与滑动平均有同样的数学表达。除有特殊情况，之后视作等同。

本章讨论频率成型滤波器的通用数学模型与程序设计。、

### 频率成型滤波器的数学模型

随机过程的频域分析法是提取信号中周期成分然后再进行积分变换在频率域上分析的一种手段。其理论基础是随机过程本身能被模拟成各种特定频率谐波叠加的系统。在这之上操作这些谐波成分的技术称为统计线性估计/滤波/平滑。

在以上假设下，将随机过程$x_t$输入以下差分系统：
$$
y_{n}=\sum_{k=0}^{M} c_{k} x_{n-k}+\sum_{j=0}^{N-1} d_{j} y_{n-j-1} \tag{1}
$$
其中$y_n$是在系统对应离散输入$x_t$的输出。通过给输入与输出设定不同的系数，就能组合出特定性质的$y_t$。当等式右边项$y$相关系数都是零时，系统只与过去时刻输入有关，具有有限项、因果性、绝对稳定性，此时称为有限递归滤波器（FIR）；当右边项$y$相关系数不为零，系统输出与过去时刻输出有关，具有条件稳定性、因果性。相比FIR系统增加了反馈环节，等同于加入了无穷多项输入，称作无穷递归滤波器（IIR）[^1]。所有的线性滤波器都是这两类滤波器的特例。

该差分系统的频率域响应：
$$
H(f)=\frac{\sum_{k=0}^{M} c_{k} e^{-2 \pi i kf }}{1-\sum_{j=0}^{N-1} d_{j} e^{-2 \pi i(j+1)f }} \tag{2}
$$
其中$f$是使用采用频率进行归一化后的频率。改变系数能相应改变频率响应函数形状，因此通过这种方式设计的滤波器称作频率成形滤波器。

特别注意在scipy.signal等软件包中，会限制$-0.5\le f\le0.5$[^3]。因为式2通过Laplace变换连续时间域得到系统传递函数，而实际上采样是离散过程，需要使用z变换到离散频率域。当离散系统采样频率$f_s\le 2f$时，Laplace变换与z变换转换中存在混叠效应。因此这些软件包中会限制一半的频率使用，如$f_s=8000Hz$，那么滤波器频率设定范围只有4000Hz。

### IIR型滤波器设计

IIR型滤波器系统比较简单：
$$
y_{n}=\sum_{k=0}^{M} c_{k} x_{n-k} \tag{3}
$$


参考文献

[^1]: 时间序列分析，汉密尔顿
[^2]: 信号与系统，奥本海姆
[^3]: scipy.signal