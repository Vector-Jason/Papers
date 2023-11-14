# -*- coding:utf-8 -*-
# @Name : ComplexMatrix
import numpy as np


# 从文件读取复数矩阵
def getComplexMatrix(M: int, N: int, file_path: str):
    Complex_Mat = np.zeros((M, N), dtype=complex)


def leftAppendRight(left: np.matrix, right: np.matrix):
    if left.shape[0] + left.shape[1] == 0:
        return np.copy(right)
    if right.shape[0] + right.shape[1] == 0:
        return np.copy(left)
    if left.shape[0] != right.shape[0]:
        return None

    return np.concatenate((left, right), axis=1)


def inverse(inputMatrix):
    if inputMatrix.shape[0] == 1:
        # tmp = np.zeros((1, 1), dtype=np.complex128)
        # tmp[0, 0] = 1 / inputMatrix[0, 0]
        # return np.matrix(tmp)
        return 1 / inputMatrix
    else:
        return np.linalg.inv(inputMatrix)