# -*- coding: utf-8 -*-

class CFileMS2:

    INDEX_ID = []
    INDEX_RT = []

    ION_MOBILITY = []

    LIST_RET_TIME = []
    ACTIVATION_CENTER = [] # isolation window中值
    WIN_SIZE = [] # 窗口大小
    LIST_PRECURSOR_ID = []  # 对应的MS1的id
    # PRECURSOR_CHARGE = []

    MATRIX_PEAK_MOZ = []  # 每一行是个list
    MATRIX_PEAK_INT = []


class CFileMS1:

    INDEX_ID = []
    INDEX_RT = []
    LIST_RET_TIME = []

    ION_MOBILITY = []

    MATRIX_PEAK_MOZ = []  # 每一行是个list
    MATRIX_PEAK_INT = []

class CSeed_FLOW2:  # 为得到Evidence准备的（碎片离子和母离子都用这个）
    MID_SCAN = []
    MOZ_FRAG = []
    INT_FRAG = []
    MID_ID = []

class CEvidence_FLOW2:

    MATRIX_PROFILE_FRAG = []
    MATRIX_MASS_DEV_FRAG = []

    LIST_RET_TIME = []   # 曲线上每个点的保留时间
    LIST_SCAN = []  # 曲线上每个点的scan

    ACTIVATION_CENTER = []
    WIN_SIZE = []
    THE_MS_INDEX = 0
    THE_MS_RT = 0
    LEFT_MID = 0
    LIST_FRAG_MOZ = []
    FRAG_MOZ = []
    FRAG_INT = []

    RT_START = 0  # 色谱曲线的起始点和终止点(存放在rt列表的位置)
    RT_END = 0
    # POINT_APEX = 0

class CEvidence_FLOW1:

    MATRIX_PROFILE_PRECURSOR = [] # 母离子profile

    LIST_RET_TIME = []   # 曲线上每个点的保留时间
    LIST_SCAN = []  # 曲线上每个点的scan
    LIST_PREC_MOZ = []

    THE_MS_INDEX = 0
    THE_MS_RT = 0 # 当前MS2对应MS1的scan号和rt
    PREC_MOZ = []
    PREC_INT = []

    RT_START = 0  # 色谱曲线的起始点和终止点(存放在rt列表的位置)
    RT_END = 0


class CDataPack:  # 这个类必须放到最后

    #  从config里面搞出来的
    LIST_PATH_MS1 = []
    LIST_PATH_MS2 = []

    LIST_PATH_ID = []

    #  需要全周期维护的
    # myID = CFileID()
    # myResult = CResult()

class CINI:

    # 这几个东东，只要爱因斯坦不被批斗，不太可能会变
    MASS_ELECTRON = 0.0005485799
    MASS_PROTON_MONO = 1.00727645224  # 1.00782503214-0.0005485799
    MASS_PROTON_ARVG = 1.0025
    # MASS_NEUTRON_AVRG = 1.0025
    MASS_NEUTRON_AVRG = 1.003

    MASS_MONO_C = 12.0
    MASS_MONO_H = 1.0078246
    MASS_MONO_N = 14.0030732
    MASS_MONO_O = 15.9949141
    MASS_MONO_S = 31.97207

    MAX_CHARGE = 4