# -*- coding:utf-8 -*-
# @Name : CompressiveReading
import os.path
import time
import joblib
import numpy as np
import pandas
import scipy.io
import matplotlib.pyplot as plt

from Signal import Signal


def leftAppendRight(left: np.matrix, right: np.matrix):
    if left.shape[0] + left.shape[1] == 0:
        return np.copy(right)
    if right.shape[0] + right.shape[1] == 0:
        return np.copy(left)
    if left.shape[0] != right.shape[0]:
        return None

    return np.concatenate((left, right), axis=1)


class CompressiveReading:
    def __init__(self):
        # 处理的信号长度
        self.N = 5000
        # 帧大小
        self.Q = 5
        # 稀疏度
        self.K = 2
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

    def generateFourierBasis(self, file_path):
        """
        生成傅里叶基
        :param file_path:
        :return:
        """
        tmp = np.matrix(np.eye(self.N, self.N), dtype=np.float64)
        phi = np.fft.fft(tmp) / np.sqrt(self.N)

        joblib.dump(phi, file_path)

    def _loadEye(self):
        """
        若不存在傅里叶基则先生成
        :return:
        """
        if self.N not in self.eyeCache:
            file_path = "./basis/{}.pkl".format(self.N)

            if not os.path.exists(file_path):
                print("Generate FourierBasis Size: {}".format(self.N))
                self.generateFourierBasis(file_path)

            self.eyeCache[self.N] = joblib.load(file_path)

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
        signal.timeIndicator = [False for _ in range(self.N)]
        # 为了消除 mod 操作的影响
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

        # 测量 / 观测矩阵
        signal.phi = np.matrix(np.zeros((self.M, self.N), dtype=np.float64))
        # 如果 timeIndicator 中为 true 这里就为 1
        for m in range(self.M):
            for n in range(m * self.Q, (m + 1) * self.Q):
                if signal.timeIndicator[n]:
                    signal.phi[m, n] = 1.0

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
        Phi: np.matrix = signal.phi
        phase: np.matrix = signal.phaseSeries

        # scipy.io.savemat("D:/Users/hoka/Projects/MATLAB/RFID/compress-sensing/DWT/phi_110.mat", dict(phi=signal.phi))
        # scipy.io.savemat("D:/Users/hoka/Projects/MATLAB/RFID/compress-sensing/DWT/phase_110.mat",
        #                  dict(phase=signal.phaseSeries))

        # 随机测量,
        s = np.dot(Phi, phase)
        # sc
        y = s.astype('complex128')
        m = 2 * self.K
        # 获取时间用于计时
        start = time.time()

        # 傅里叶基
        fftEyeNN = np.matrix(self._loadEye())

        # 稀疏矩阵
        Psi = fftEyeNN * (1 / np.sqrt(self.N))
        # Psi 的共轭转置, trans
        Psi_H = Psi.H

        # 传感 / 感知矩阵, T
        # Θ = Φ Ψ^(-1)
        # 恢复矩阵(测量矩阵*正交反变换矩阵)
        Theta = np.matrix(np.zeros((self.M, self.N), dtype=np.complex128))
        # Θ = Φ Ψ^(-1)
        for i in range(self.M):
            for j in range(self.N):
                # 前为 实部, 后为 虚部
                result = [0.0, 0.0]
                for k in range(i*self.Q, (i+1)*self.Q):
                    if Phi[i, k] == 1:
                        tmp = Psi_H[k, j]
                        result[0] += tmp.real
                        result[1] += tmp.imag
                Theta[i, j] = result[0] + 1j * result[1]

        # 待重构的谱域(变换域)向量, hat_y
        S = np.matrix(np.zeros((1, self.N), dtype=np.complex128))
        # 增量矩阵(初始值为空矩阵), Aug_t
        Theta_new = np.matrix(np.empty((0, 0), dtype=np.complex128))
        # Aug_y
        Lambda = None
        # 残差值, 初始为 s, r_n
        r = s.astype('complex128')
        # 存放贡献值, 计算 感知矩阵Theta 对 残差r 的贡献, product
        omega = [0.0 for _ in range(self.N)]
        # 存储 m 次迭代的最大值列索引
        pos_array = [0 for _ in range(m)]

        # 零列, 用于填充
        zeroColumn = np.matrix(np.zeros((self.M, 1), dtype=np.complex128))

        for times in range(m):
            for col in range(self.N):
                # 计算 感知矩阵Theta 对 残差r 的贡献
                temp = np.dot(Theta[:, col].H, r)
                omega[col] = np.sqrt(temp[0, 0].real ** 2 + temp[0, 0].imag ** 2)

            # 选择 omega 中最大贡献值
            # pos 为傅里叶域最大分量的频率下标, 最大投影系数对应的位置
            pos = int(np.argmax(omega))
            print("pos / {}".format(self.N / 1000), pos / (self.N / 1000))

            # 矩阵扩充, 从 1 到 10, 每次都加上最大值所在列
            Theta_new = np.matrix(leftAppendRight(Theta_new, Theta[:, pos]))
            # 将选中的列去掉，为了简单把它置零
            Theta[:, pos] = zeroColumn
            # 计算 Theta_new 对 残差r 的贡献, Aug_y
            # 广义逆矩阵（伪逆）(A^T*A)^(-1)*A^T，最小二乘法得到的解和伪逆矩阵是等价的
            middle_res = np.dot(np.linalg.inv(np.dot(Theta_new.H, Theta_new)), Theta_new.H)
            # 计算实际的贡献, 即最小化 min‖Theta_new Lambda - y‖_2
            # 最小二乘,使残差最小
            Lambda = np.dot(middle_res, y)
            # 计算剩余残差
            r = y - np.dot(Theta_new, Lambda)
            # 记录这轮的最大系数位置
            # 纪录最大投影系数的位置
            pos_array[times] = pos

        for i in range(m):
            # 重构向量
            # 重构的谱域向量
            S[0, pos_array[i]] = Lambda[i, 0]

        # 已知 y、Φ、Ψ, 求出 S(即代码中的 hat_y), 再用 x = ΨS 求出高维信号 x 就为最后恢复的信号
        # 做逆傅里叶变换重构得到时域信号
        x = np.dot(Psi.H, S.T)
        signal.recoveredSeries = x.real
        print("CompressiveReading time: ", time.time() - start)

        return signal


if __name__ == '__main__':
    cr = CompressiveReading()
    df = pandas.read_csv("./data/exp-340hz-90s.csv")
    # df = pandas.read_csv("./data/exp-340hz-60s.csv")
    i = 1
    phaseList = np.array(df["Phase"], dtype=np.float64)[i * cr.N:(i + 1) * cr.N]
    timeList = np.array(df["FirstReadTime"])[i * cr.N:(i + 1) * cr.N]
    # phaseList = np.array(df["Phase"], dtype=np.float64)[i * 1000:(i + 1) * 1000]
    # timeList = np.array(df["FirstReadTime"])[i * 1000:(i + 1) * 1000]
    phaseList = phaseList / 4096.0 * 2 * np.pi

    # TODO 放大微振动
    # 非直流分量
    # B = phaseList - np.mean(phaseList)
    # phaseList += 150 * B
    # phaseList = phaseList % (2 * np.pi)

    phaseList = phaseList - np.mean(phaseList)

    print(len(timeList), len(phaseList))

    signal = cr.recover(timeList, phaseList)
    output = np.squeeze(np.array(signal.recoveredSeries[:, 0]))
    output_phase = np.squeeze(np.array(signal.phaseSeries[:, 0]))
    # joblib.dump(output, "D:/Users/hoka/Projects/Python/RFID/tagbeat/data/138hz-90s/recover_1.pkl")
    # joblib.dump(output_phase, "D:/Users/hoka/Projects/Python/RFID/tagbeat/data/138hz-90s/phase_1.pkl")
    print(output.shape)

    sampling_rate = 1000  # 由文件算出
    fft_size = cr.N

    xs = output[:fft_size]
    # xs = xs - np.mean(xs)
    # scipy.io.savemat('D:/Users/hoka/Projects/MATLAB/RFID/compress-sensing/powerSpectrum/output_340.mat', dict(x=xs))
    xf = np.fft.rfft(xs) / fft_size
    freqs = np.linspace(0, int(sampling_rate / 2), int(fft_size / 2) + 1)
    xfp = np.abs(xf)
    # 去除 0-1 hz 内的分量
    xfp[0: 6] = 0
    print(freqs[np.argmax(xfp)])

    plt.figure(figsize=(16, 8))
    plt.subplot(211)
    plt.plot(xs)
    plt.xlabel("data")
    plt.subplot(212)
    plt.plot(freqs, xfp)
    plt.xlabel(u"Hz")
    plt.subplots_adjust(hspace=0.4)
    plt.show()
