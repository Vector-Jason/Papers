# -*- coding:utf-8 -*-
# @Name : test
# coding:utf-8
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DCT基作为稀疏基，重建算法为OMP算法 ，图像按列进行处理
# 参考文献: 任晓馨. 压缩感知贪婪匹配追踪类重建算法研究[D].
# 北京交通大学, 2012.
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 导入所需的第三方库文件
import numpy as np
import math
from PIL import Image
import matplotlib.pyplot as plt


def cs_omp(y, Phi, K):
    residual = y  # 初始化残差
    (M, N) = Phi.shape
    index = np.zeros(N, dtype=int)
    for i in range(N):  # 第i列被选中就是1，未选中就是-1
        index[i] = -1
    result = np.zeros((N, 1))
    for j in range(K):  # 迭代次数
        # 计算各原子对 y 的贡献值
        product = np.fabs(np.dot(Phi.T, residual))
        pos = np.argmax(product)  # 最大投影系数对应的位置
        index[pos] = 1  # 对应的位置取1
        my = np.linalg.pinv(Phi[:, index >= 0])  # 广义逆矩阵（伪逆）(A^T*A)^(-1)*A^T，最小二乘法得到的解和伪逆矩阵是等价的。
        a = np.dot(my, y)  # 最小二乘
        # 剩余的残差
        residual = y - np.dot(Phi[:, index >= 0], a)
    result[index >= 0] = a
    Candidate = np.where(index >= 0)  # 返回所有选中的列
    return result, Candidate


if __name__ == '__main__':
    # 读取图像，并变成numpy类型的 array
    im = Image.open('lena.jpg').convert('L')  # 模式“RGB”转换为其他不同模式
    # 图片大小256*256，为灰度图像，每个像素用8个bit表示，0表示黑，255表示白，其他数字表示不同的灰度。
    # 转换公式：L = R * 299/1000 + G * 587/1000+ B * 114/1000。
    plt.imshow(im, cmap='gray')  # cmap参数: 为调整显示颜色  gray为黑白色，加_r取反为白黑色
    N = 256
    M = 128
    K = 50
    im = np.array(im)
    # 观测矩阵
    Phi = np.random.randn(M, N) / np.sqrt(M)
    # 随机测量
    img_cs_1d = np.dot(Phi, im)  # （M*N）

    # 生成稀疏基DCT矩阵
    # TODO DCT基矩阵
    mat_dct_1d = np.zeros((N, N))
    v = range(N)
    for k in range(0, N):
        dct_1d = np.cos(np.dot(v, k * math.pi / N))
        if k > 0:
            dct_1d = dct_1d - np.mean(dct_1d)
        mat_dct_1d[:, k] = dct_1d / np.linalg.norm(dct_1d)
    print(mat_dct_1d)

    # 重建
    sparse_rec_1d = np.zeros((N, N))  # 初始化稀疏系数矩阵
    # M x N * N x N = M x N
    Theta_1d = np.dot(Phi, mat_dct_1d)  # 测量矩阵乘上稀疏基矩阵

    for i in range(N):
        if i % 32 == 0:
            print('正在重建第', i, '列。')
        y = np.reshape(img_cs_1d[:, i], (M, 1))  # 将图像的第I列进行重构
        column_rec, Candidate = cs_omp(y, Theta_1d, K)  # 利用OMP算法计算稀疏系数，column_rec为重构的稀疏系数的值，不是x的值，x=稀疏基与稀疏系数的值
        # 做内积
        x_pre = np.reshape(column_rec, (N))  # 常用于矩阵规格变换，将矩阵转换为特定的行和列的矩阵
        sparse_rec_1d[:, i] = x_pre
    img_rec = np.dot(mat_dct_1d, sparse_rec_1d)  # 稀疏系数乘上基矩阵，对图像进行重构
    # 显示重建后的图片
    img_pre = Image.fromarray(img_rec)
    plt.imshow(img_pre, cmap='gray')
    plt.show()
    error = np.linalg.norm(img_rec - im) / np.linalg.norm(im)  # 采样率为M/N=50%,说明在余弦变换下或许不满足稀疏度假设，也就是在此稀疏基下不稀疏
    error
