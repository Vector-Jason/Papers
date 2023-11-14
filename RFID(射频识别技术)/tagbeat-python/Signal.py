# -*- coding:utf-8 -*-
# @Name : Signal
class Signal:
    def __init__(self):
        self.CODE_NO_ERROR = 0
        self.CODE_NO_VIBRATION = 1
        self.code = None

        self.timeIndicator = None
        # 相位序列
        self.phaseSeries = None
        # IV.C 中的测量矩阵
        self.phi = None
        # 恢复序列
        self.recoveredSeries = None
        # 基频
        self.frequency = None
        # 原始时间序列, 注意其长度和上面的不同
        self.timeList = None
        # 原始相位序列, 注意其长度和上面的不同
        self.originPhaseList = None
        # timeList 中用到的最后一个时间索引
        self.lastTimeIndex = 0


