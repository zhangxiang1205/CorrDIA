# -*- coding: utf-8 -*-
from MSData import CFileMS1, CFileMS2, CSeed_FLOW2, CEvidence_FLOW2, CEvidence_FLOW2, CEvidence_FLOW1
from MSOperator import opGetStartAndEndForProfile
import pandas
import json
import numpy as np

from MSSystem import VALUE_MAX_SCAN
from MSTool import toolFindNeighborFromSortedList1, toolFindIndexFromSortedList1


class CFunctionEvidence():

    # def __captainGetCycleWinNum2(self, inputDataMS2, MidIndex):
    #
    #     centerMoz = inputDataMS2.LIST_ACTIVATION_CENTER[MidIndex]
    #     ion_mob = inputDataMS2.ION_MOBILITY[MidIndex]
    #     try:
    #         num = 1
    #         while centerMoz != inputDataMS2.LIST_ACTIVATION_CENTER[MidIndex + num] or \
    #                 ion_mob != inputDataMS2.ION_MOBILITY[MidIndex + num]:
    #             num = num + 1
    #     except IndexError:  # 超出激活中心列表的边界，从左边计算
    #         num = 1
    #         while centerMoz != inputDataMS2.LIST_ACTIVATION_CENTER[MidIndex - num]:
    #             num = num + 1
    #
    #     return num
    def __captainGetCycleWinNum(self, inputDataMS2: CFileMS2, MidIndex):

        centerMoz = inputDataMS2.ACTIVATION_CENTER[MidIndex]
        try:
            num = 1
            while centerMoz != inputDataMS2.ACTIVATION_CENTER[MidIndex + num]:
                num = num + 1
        except IndexError:  # 超出激活中心列表的边界，从左边计算
            num = 1
            while centerMoz != inputDataMS2.ACTIVATION_CENTER[MidIndex - num]:
                num = num + 1

        return num

    # def __captainGetCycleWinNum1(self, inputDataMS1, MidIndex):
    #
    #     ion_mob = inputDataMS1.ION_MOBILITY[MidIndex]
    #     try:
    #         num = 1
    #         while ion_mob != inputDataMS1.ION_MOBILITY[MidIndex + num]:
    #             num = num + 1
    #     except IndexError:  # 超出激活中心列表的边界，从左边计算
    #         num = 1
    #         while ion_mob != inputDataMS1.ION_MOBILITY[MidIndex - num]:
    #             num = num + 1
    #
    #     return num

    def __captainGetIndexListFromMS2(self, inputDataMS2, Midindex, inputWinRT):

        cycle =  self.__captainGetCycleWinNum(inputDataMS2, Midindex) # 遍历一个cycle有多少个ms2(我目前的数据是64)
        # cycle = 22 # 遍历一个cycle有多少个ms2(我目前的数据是64，cycle=65)
        MidRT = inputDataMS2.LIST_RET_TIME[Midindex] # 找RT(取左右2min范围内)

        RT_left = MidRT
        tmp_index = Midindex
        left_list_ms2 = []
        while (RT_left >= MidRT - inputWinRT) & (tmp_index >= cycle):
            tmp_index -= cycle
            if tmp_index >1:
                RT_left = inputDataMS2.LIST_RET_TIME[tmp_index]
                left_list_ms2.append(tmp_index)
        left_list_ms2.reverse()

        RT_right = MidRT
        tmp_index = Midindex
        right_list_ms2 = []
        max_time = inputDataMS2.INDEX_RT[-1]
        idx = inputDataMS2.LIST_RET_TIME.index(max_time)
        while (RT_right <= MidRT + inputWinRT) & (tmp_index + cycle <= idx):
            tmp_index += cycle
            if tmp_index > 1:
                RT_right = inputDataMS2.LIST_RET_TIME[tmp_index]
                right_list_ms2.append(tmp_index)

        index_list_ms2 = left_list_ms2 + [Midindex] + right_list_ms2
        iMid = len(left_list_ms2)

        # length = len(index_list_ms2) # MS2保留时间范围内的list长度
        return index_list_ms2, iMid, cycle, left_list_ms2, right_list_ms2

    def __captainGetIndexListFromMS1(self, inputDataMS1, Midindex, inputWinRT, left_list_ms2, right_list_ms2 ):

        cycle = 22 # 遍历一个cycle有多少个ms1(我目前的数据是64)
        # cycle = self.__captainGetCycleWinNum1(inputDataMS1, Midindex)
        # MidRT = inputDataMS1.LIST_RET_TIME[Midindex] # 找RT(取左右2min范围内)
        #
        # RT_left = MidRT
        # tmp_index = Midindex
        # left_list_ms1 = []
        # while (RT_left >= MidRT - inputWinRT) & (tmp_index >= cycle):
        #     tmp_index -= cycle
        #     if tmp_index > 0:
        #         RT_left = inputDataMS1.LIST_RET_TIME[tmp_index]
        #         left_list_ms1.append(tmp_index)
        # left_list_ms1.reverse()
        #
        # RT_right = MidRT
        # tmp_index = Midindex
        # right_list_ms1 = []
        # while (RT_right <= MidRT + inputWinRT) & (tmp_index + cycle <len(inputDataMS1.LIST_RET_TIME)):
        #     tmp_index += cycle
        #     RT_right = inputDataMS1.LIST_RET_TIME[tmp_index]
        #     right_list_ms1.append(tmp_index)
        #
        # index_list_ms1 = left_list_ms1 + [Midindex] + right_list_ms1
        # iMid = len(left_list_ms1)
        left_list_ms1 = []
        right_list_ms1 = []
        if len(left_list_ms2) == 0: left_list_ms1 = []
        if(len(right_list_ms2) == 0): right_list_ms1 = []
        for i, data in enumerate(left_list_ms2):
            tmp = left_list_ms2[i] - (left_list_ms2[i] - 1) % cycle
            left_list_ms1.append(tmp)
        for i, data in enumerate(right_list_ms2):
            tmp = right_list_ms2[i] - (right_list_ms2[i] - 1) % cycle
            right_list_ms1.append(tmp)
        # new_Midindex = Midindex - (Midindex - 1) % cycle
        index_list_ms1 = left_list_ms1 + [Midindex] + right_list_ms1
        iMid = len(left_list_ms1)
        return index_list_ms1, iMid

    def __captainfillEvidence1(self, the_index_ms1, index_list_ms1, mid, inputdataMS1: CFileMS1,inputSeed: CSeed_FLOW2, flagGetStartAndEnd: True):

        listScan = index_list_ms1 # 保留时间在2min之内的MS1索引列表
        print(listScan)
        nScan = len(listScan)
        accuracy = 0.02
        listPrec_MOZ = inputSeed.MOZ_FRAG
        # for moz in list_MOZ:
        #     if int(float(moz)) <= 1025 and int(float(moz)) >= 1000:
        #         listPrec_MOZ.append(moz)

        outputEvidence = CEvidence_FLOW1()
        # outputEvidence.LIST_PREC_MOZ = [[]] * len(listPrec_MOZ)
        # outputEvidence.LIST_PREC_MOZ = np.zeros(shape=[len(listPrec_MOZ), nScan])
        outputEvidence.PREC_MOZ = inputSeed.MOZ_FRAG
        outputEvidence.PREC_INT = inputSeed.INT_FRAG
        outputEvidence.MATRIX_PROFILE_PRECURSOR = np.zeros(shape=[len(listPrec_MOZ), nScan])
        outputEvidence.LIST_PREC_MOZ = np.zeros(shape=[len(listPrec_MOZ), nScan])
        outputEvidence.LIST_SCAN = [[]] * nScan
        outputEvidence.LIST_RET_TIME = [[]] * nScan
        outputEvidence.POINT_APEX = mid
        outputEvidence.THE_MS_INDEX = 0
        outputEvidence.THE_MS_RT = 0

        # for i, data in enumerate(listPrec_MOZ):
        #     outputEvidence.LIST_PREC_MOZ[i] = data

        for i, iScan in enumerate(listScan):
            # print(i, iScan)
            tmp_MOZ_list = inputdataMS1.MATRIX_PEAK_MOZ[iScan]
            tmp_INT_list = inputdataMS1.MATRIX_PEAK_INT[iScan]

            outputEvidence.LIST_SCAN[i] = iScan
            outputEvidence.LIST_RET_TIME[i] = inputdataMS1.LIST_RET_TIME[iScan]
            outputEvidence.THE_MS_INDEX = the_index_ms1
            outputEvidence.THE_MS_RT = inputdataMS1.LIST_RET_TIME[the_index_ms1]

            for j, fragment_moz in enumerate(listPrec_MOZ):
                indexINT = toolFindNeighborFromSortedList1(tmp_MOZ_list, fragment_moz)

                mass_dev = abs(float(tmp_MOZ_list[indexINT]) - float(fragment_moz))
                if mass_dev < accuracy:
                    outputEvidence.LIST_PREC_MOZ[j, i] = tmp_MOZ_list[indexINT]
                    outputEvidence.MATRIX_PROFILE_PRECURSOR[j, i] = tmp_INT_list[indexINT]
                else:
                    outputEvidence.LIST_PREC_MOZ[j, i] = 0.0
                    outputEvidence.MATRIX_PROFILE_PRECURSOR[j, i] = 0.0

        frag_rank_list = np.argsort(-1 * outputEvidence.MATRIX_PROFILE_PRECURSOR[:, mid])
        outputEvidence.MATRIX_PROFILE_PRECURSOR = outputEvidence.MATRIX_PROFILE_PRECURSOR[frag_rank_list]
        outputEvidence.LIST_PREC_MOZ = outputEvidence.LIST_PREC_MOZ[frag_rank_list]
        PREC_MOZ = []
        PREC_INT = []
        for i in frag_rank_list:
            PREC_MOZ.append(outputEvidence.PREC_MOZ[i])
            PREC_INT.append(outputEvidence.PREC_INT[i])
        outputEvidence.PREC_MOZ = PREC_MOZ
        outputEvidence.PREC_INT = PREC_INT
        # 得到色谱曲线的起始点和终止点
        rt_start_list = []
        rt_end_list = []
        if len(listPrec_MOZ) > 0:
            for i in range(outputEvidence.MATRIX_PROFILE_PRECURSOR.shape[0]):
                i_left, i_right = opGetStartAndEndForProfile(outputEvidence.MATRIX_PROFILE_PRECURSOR[i], mid, 0.1, 1)
                rt_start_list.append(i_left)
                rt_end_list.append(i_right)
        rtStartPosition = np.argmax(np.bincount(rt_start_list))
        rtEndPosition = np.argmax(np.bincount(rt_end_list))
        outputEvidence.RT_START = rtStartPosition
        outputEvidence.RT_END = rtEndPosition

        return outputEvidence

    def __captainfillEvidence2(self, index_list_ms2, mid, inputdataMS2: CFileMS2,inputSeed: CSeed_FLOW2, the_index_ms2, flagGetStartAndEnd: True):
        listScan = index_list_ms2 # 保留时间在范围之内的MS2索引列表
        print(listScan)
        nScan = len(listScan)
        accuracy = 0.02
        listFragment_MOZ = inputSeed.MOZ_FRAG # 这张MS2所有的谱峰质荷比信息
        # init

        outputEvidence = CEvidence_FLOW2()
        outputEvidence.FRAG_MOZ = inputSeed.MOZ_FRAG
        outputEvidence.FRAG_INT = inputSeed.INT_FRAG
        outputEvidence.MATRIX_PROFILE_FRAG = np.zeros(shape=[len(listFragment_MOZ), nScan])
        outputEvidence.LIST_FRAG_MOZ = np.zeros(shape=[len(listFragment_MOZ), nScan])
        outputEvidence.LIST_SCAN = [[]] * nScan
        outputEvidence.LIST_RET_TIME = [[]] * nScan
        outputEvidence.THE_MS_INDEX = the_index_ms2
        outputEvidence.LEFT_MID = mid

        outputEvidence.THE_MS_RT = inputdataMS2.LIST_RET_TIME[the_index_ms2]
        outputEvidence.ACTIVATION_CENTER = inputdataMS2.ACTIVATION_CENTER[the_index_ms2]
        outputEvidence.WIN_SIZE = inputdataMS2.WIN_SIZE[the_index_ms2]

        for i, iScan in enumerate(listScan):

            tmp_MOZ_list = inputdataMS2.MATRIX_PEAK_MOZ[iScan]
            tmp_INT_list = inputdataMS2.MATRIX_PEAK_INT[iScan]

            outputEvidence.LIST_SCAN[i] = iScan
            outputEvidence.LIST_RET_TIME[i] = inputdataMS2.LIST_RET_TIME[iScan]

            for j, fragment_moz in enumerate(listFragment_MOZ):
                indexINT = toolFindNeighborFromSortedList1(tmp_MOZ_list, fragment_moz)

                mass_dev = abs(float(tmp_MOZ_list[indexINT]) - fragment_moz)
                if mass_dev < accuracy:
                    outputEvidence.LIST_FRAG_MOZ[j, i] = tmp_MOZ_list[indexINT]
                    outputEvidence.MATRIX_PROFILE_FRAG[j, i] = tmp_INT_list[indexINT]
                else:
                    # 误差超过20ppm
                    outputEvidence.LIST_FRAG_MOZ[j, i] = 0.0
                    outputEvidence.MATRIX_PROFILE_FRAG[j, i] = 0.0

        # 按照mid这一列强度降序排列的顺序更新整个矩阵
        frag_rank_list = np.argsort(-1 * outputEvidence.MATRIX_PROFILE_FRAG[:, mid])
        outputEvidence.MATRIX_PROFILE_FRAG = outputEvidence.MATRIX_PROFILE_FRAG[frag_rank_list]
        outputEvidence.LIST_FRAG_MOZ = outputEvidence.LIST_FRAG_MOZ[frag_rank_list]
        FRAG_MOZ = []
        FRAG_INT = []
        for i in frag_rank_list:
            FRAG_MOZ.append(outputEvidence.FRAG_MOZ[i])
            FRAG_INT.append(outputEvidence.FRAG_INT[i])
        outputEvidence.FRAG_MOZ = FRAG_MOZ
        outputEvidence.FRAG_INT = FRAG_INT
        # 得到色谱曲线的起始点和终止点
        rt_start_list = []
        rt_end_list = []
        if len(listFragment_MOZ) > 0:
            for i in range(outputEvidence.MATRIX_PROFILE_FRAG.shape[0]):
                i_left, i_right = opGetStartAndEndForProfile(outputEvidence.MATRIX_PROFILE_FRAG[i], mid, 0.1, 1)
                rt_start_list.append(i_left)
                rt_end_list.append(i_right)
        rtStartPosition = np.argmax(np.bincount(rt_start_list))
        rtEndPosition = np.argmax(np.bincount(rt_end_list))
        outputEvidence.RT_START = rtStartPosition
        outputEvidence.RT_END = rtEndPosition

        return outputEvidence

    def fillEvidence1(self, inputdataMS1: CFileMS1, iID, inputSeed: CSeed_FLOW2, left_list_ms2, right_list_ms2 ):

        winRT = 0.5 # 变量 RT氛围1min之内的窗口

        iMid = iID  # 当前MS1的scan号

        # 得到保留时间窗口范围内的索引列表
        index_list_ms1, mid = self.__captainGetIndexListFromMS1(inputdataMS1, iMid, winRT, left_list_ms2, right_list_ms2 )

        return self.__captainfillEvidence1(iID, index_list_ms1, mid, inputdataMS1, inputSeed, True)

    def fillEvidence2(self, inputdataMS2: CFileMS2, iID, inputSeed: CSeed_FLOW2):

        winRT = 0.5 # 变量 RT氛围1min之内的窗口

        iMid = iID
        # 假设当前只画id为1的1000-1025窗口内的碎片离子色谱曲线,在dataMS2中为第一个

        # 得到保留时间窗口范围内的索引列表
        index_list_ms2, mid, cycle, left_list_ms2, right_list_ms2 = self.__captainGetIndexListFromMS2(inputdataMS2, iMid, winRT)

        return self.__captainfillEvidence2(index_list_ms2, mid, inputdataMS2, inputSeed, iID, True), cycle, index_list_ms2, left_list_ms2, right_list_ms2