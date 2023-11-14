# -*- coding:utf-8 -*-
# @Name : tagbeatFunc
import math
import os.path
import time
import joblib
import numpy as np
import pandas

from Signal import Signal
from ComplexMatrix import leftAppendRight, inverse


class CompressiveReading:
    def __init__(self):
        # 处理的信号长度
        self.N = 5000
        # 帧大小
        self.Q = 5
        # 稀疏度
        self.K = 5
        # 行数
        self.M = int(np.ceil(self.N / self.Q))

        # 傅里叶基
        self.eyeCache = {}

    def changeSampleNumber(self, N):
        self.N = N
        self.M = int(np.ceil(self.N / self.Q))
        print("Change parameter of N to ", N)

    def changeFrameSize(self, Q):
        self.Q = Q
        self.M = int(np.ceil(self.N / self.Q))
        print("Change parameter of Q to", Q)

    def changeSparsity(self, K):
        self.K = K
        print("Change parameter of K to", K)

    def _loadEye(self):
        if self.N not in self.eyeCache:
            file_path = "./basis/{}.pkl".format(self.N)

            if not os.path.exists(file_path):
                raise Exception('[ERROR]:Cannot find the file of Fourier bases.')

            # self.eyeCache[self.N] = np.fromfile(file_path)
            self.eyeCache[self.N] = joblib.load(file_path)

            # print(self.eyeCache[self.N])
            # print(self.eyeCache[self.N][0][0])
            # print(self.eyeCache[self.N].shape)
            print(self.eyeCache[self.N].dtype)

        return self.eyeCache[self.N]

    def _preprocess(self, time: list, phase: list):

        print("N[{}], Q[{}], K[{}]".format(self.N, self.Q, self.K))

        # 确保样本长度大于 N
        base = time[0]
        for i in range(len(time)):
            # 缩小 1000 倍, 微秒转毫秒
            time[i] = np.floor(time[i] - base) / 1000.0

        signal = Signal()
        # 采样矩阵
        signal.timeIndicator = [False for i in range(self.N)]
        # N x 1 的取 sin 后的相位序列
        # > 为了消除 mod 操作的影响
        signal.phaseSeries = np.matrix(np.zeros((self.N, 1), dtype=np.float64))

        # 取 N 个值, 相当于 5 秒内的数据
        for i in range(len(time)):
            index = int(time[i])
            if index < self.N:
                signal.timeIndicator[index] = True
                signal.phaseSeries[index, 0] = np.sin(phase[i])
                # 记录下最后一个时间索引
                signal.lastTimeIndex = i
            else:
                break

        # IV.C 中的测量矩阵
        signal.phi = np.matrix(np.zeros((self.M, self.N), dtype=np.float64))
        # 如果 timeIndicator 中为 true 这里就为 1
        for m in range(self.M):
            for n in range(m * self.Q, (m + 1) * self.Q):
                if signal.timeIndicator[n]:
                    signal.phi[m, n] = 1.0
        # TODO 打印下 phi 是否如下形式
        """
        1, 1, 0, 0, 0, 0
        0, 0, 1, 1, 0, 0,
        0, 0, 0, 0, 1, 1
        """

        return signal

    def recover(self, time, phase):
        # 输入的 time, phase 长度为 1000, 输出长度 5000 的 phase 和 recoveredPhase
        tmp_time = np.copy(time)
        tmp_phase = np.copy(phase)

        signal = self._preprocess(time, phase)

        signal.timeList = tmp_time
        signal.originPhaseList = tmp_phase

        return self._seekFrequency(self._recover(signal))

    def _seekFrequency(self, signal):
        return signal

    def _recover(self, signal):
        phi: np.matrix = signal.phi
        phase: np.matrix = signal.phaseSeries

        # 实矩阵: M x 1 = (M x N) .* (N x 1)
        # s = phi * phase = 论文中的 y ?
        # 测量值
        s: np.matrix = np.dot(phi, phase)
        # 复数矩阵: s 为实部
        sc = s.astype('complex128')

        # 10
        m = 2 * self.K
        # 获取当前时间, 用于计时
        start = time.time()

        # 傅里叶基: 复数矩阵
        # N * N
        fftEyeNN = np.matrix(self._loadEye())

        # 除以 sqrt(5000)
        # N * N
        Psi = fftEyeNN * (1 / np.sqrt(self.N))

        # Psi 的共轭转置
        # N x N
        trans = Psi.H

        # M x N : 1000 x 5000
        # 随机矩阵与稀疏表示基的结合
        T = np.matrix(np.zeros((self.M, self.N), dtype=np.complex128))
        # T = phi * psi^(-1)
        # 对应 Θ = ΦΨ, 称为传感矩阵, y = ΘS
        for i in range(self.M):
            for j in range(self.N):
                # 前为 实部, 后为 虚部
                result = [0.0, 0.0]
                for k in range(i*self.Q, (i+1)*self.Q):
                    if phi[i, k] == 1:
                        tmp = trans[k, j]
                        result[0] += tmp.real
                        result[1] += tmp.imag
                T[i, j] = result[0] + 1j * result[1]

        # 1 x N, 是否是 稀疏系数 S?
        hat_y = np.matrix(np.zeros((1, self.N), dtype=np.complex128))
        # 增量矩阵 (初始为空)
        Aug_t = np.matrix(np.empty((0, 0), dtype=np.complex128))
        Aug_y = None
        # 残差值, 初始残差为 s
        r_n = s.astype('complex128')
        # temp = None
        # 1 x N = 1 x 5000
        product = [0.0 for _ in range(self.N)]
        # m = 2 x K = 10
        # 存储 pos: product 中最大值所在 T 列
        pos_array = [0 for _ in range(m)]
        # M x 1 = 1000 x 1
        zeroColumn = np.matrix(np.zeros((self.M, 1), dtype=np.complex128))
        # zeroColumn = [0.0 for _ in range(self.M)]

        for times in range(m):
            # TODO: 该部分可优化, 可以只运算一次
            for col in range(self.N):
                # 1 x 1 = (1 x M) x (M x 1)
                temp = np.dot(T[:, col].H, r_n)
                product[col] = np.sqrt(temp[0, 0].real ** 2 + temp[0, 0].imag ** 2)

            # 最大值的索引
            # 这个 pos 就是傅里叶域最大分量的频率下标
            # 最大投影系数对应的位置
            pos = int(np.argmax(product))
            print("pos / 3", pos / 3)

            # Aug_t = np.concatenate((Aug_t, T[:, pos]), axis=1)
            # Aug_t 逐渐累积, 从 1 到 10, 每次都加上最大值所在列
            # 矩阵扩充
            Aug_t = np.matrix(leftAppendRight(Aug_t, T[:, pos]))
            # 选中的列置零
            T[:, pos] = zeroColumn
            # (1 x 1) = (1 x 1000) x (1000 x 1)
            # (2 x 2) = (2 x 1000) x (1000 x 2)
            # :
            # :
            # (10 x 10) = (10 x 1000) x (1000 x 10)
            middle_res = np.dot(Aug_t.H, Aug_t)
            # (1 x 1) = [(1 x 1) x (1 x 1000)] x (1000 x 1)
            # (2 x 1) = [(2 x 2) x (2 x 1000)] x (1000 x 1)
            # :
            # :
            # (10 x 1) = [(10 x 10) x (10 x 1000)] x (1000 x 1)
            # 最小二乘
            Aug_y = np.dot(np.dot(np.linalg.inv(middle_res), Aug_t.H), sc)
            # (1000 x 1) = (1000 x 1) - [(1000 x 1) x (1 x 1)]
            # (1000 x 1) = (1000 x 1) - [(1000 x 2) x (2 x 1)]
            # :
            # :
            # (1000 x 1) = (1000 x 1) - [(1000 x 10) x (10 x 1)]
            # 计算残差
            r_n = sc - np.dot(Aug_t, Aug_y)
            # 记录最大投影系数的位置
            pos_array[times] = pos

        for i in range(m):
            # 1 x 5000 的 [0, pos] = (10 x 1) [0, 0], [1, 0] ... [9, 0]
            print("Aug_y", Aug_y[i, 0])
            # 重构的向量
            hat_y[0, pos_array[i]] = Aug_y[i, 0]

        # (5000 x 1) = (5000 x 5000) x (5000 x 1)
        # x(未知) = \Psi(稀疏矩阵) x S(稀疏系数)
        # y = ΦΨS
        # 已知 y、Φ、Ψ, 求出 S(即代码中的 hat_y), 再用 x = ΨS 求出高维信号 x 就为最后恢复的信号
        # 待重构的谱域 (变换域) 向量
        hat_x = np.dot(Psi.H, hat_y.T)
        # 最后的信号是 hat_x
        signal.recoveredSeries = hat_x.real
        print("CompressiveReading time: ", time.time() - start)

        return signal


if __name__ == '__main__':
    cr = CompressiveReading()
    # df = pandas.read_csv("D:/Users/hoka/Projects/Java/RFID/tagbeat/benchmark/test/exp-15hz-30s.csv")
    df = pandas.read_csv("D:/Users/hoka/Projects/Java/RFID/tagbeat/benchmark/test/exp-275hz-90s.csv")
    i = 2
    phaseList = np.array(df["Phase"], dtype=np.float64)[i * 5000:(i + 1) * 5000]
    timeList = np.array(df["FirstReadTime"])[i * 5000:(i + 1) * 5000]
    phaseList = phaseList / 4096.0 * 2 * np.pi

    # TODO 放大微振动
    # 非直流分量
    B = phaseList - np.mean(phaseList)
    phaseList += 0.2 * B
    phaseList = phaseList % (2 * np.pi)

    print(len(timeList), len(phaseList))

    signal = cr.recover(timeList, phaseList)
    output = np.squeeze(np.array(signal.recoveredSeries[:, 0]))
    joblib.dump(output, "D:/Users/hoka/Projects/Python/RFID/FFA_riptide/data/275hz-90s/recover_2.pkl")
    print(output.shape)

