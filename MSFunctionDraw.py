#  重构色谱曲线
import multiprocessing
import time

import pandas
import matplotlib.pyplot as plt

import MSFunctionSimilar
from MS1Preprocess import GetCharge
from MSFunctionEvidence import CFunctionEvidence
from MSFunctionSimilar import CFunctionSimilar
from MSTool import CFunctionParseMS1, CFunctionParseMS2, toolFindNeighborFromSortedList1
from MSData import CFileMS2, CFileMS1, CSeed_FLOW2
import numpy as np
from MSSystem import VALUE_ILLEGAL, VALUE_MAX_SCAN


class CFunctionDrawDIAXIC:

    def __captainFillSeed(self, iID, dataMS2: CFileMS2):

        # 初始化seed
        myseed = CSeed_FLOW2()
        myseed.MOZ_FRAG = dataMS2.MATRIX_PEAK_MOZ[iID]
        myseed.INT_FRAG = dataMS2.MATRIX_PEAK_INT[iID]

        return myseed

    def __captainGetEvidence(self, LIST_PATH_MS1, LIST_PATH_MS2, inputListEvidence1, inputListSeed1, inputListEvidence,
                             inputListSeed, inputListIndex):
        '''
        :param inputListEvidence: 输入空的evidence列表用于填充，提取强度点，获取色谱曲线矩阵用于画图
        :param inputListSeed: 输入空的seed列表用于填充，得到的seed用于对应evidence的获取
        :param inputListIndex: 输入索引列表，去数据结构中寻址得到rt等信息
        :return:
        '''
        # 全部MS2做拆谱，即inputListIndex里面为所有MS2的scan号
        global cycle
        # the_index_start = LIST_PATH_MS2.INDEX_ID.index(14170)
        # the_index_end = LIST_PATH_MS2.INDEX_ID.index(14190)
        # inputListIndex = LIST_PATH_MS2.INDEX_ID[the_index_start: the_index_end+1]

        the_index = LIST_PATH_MS2.INDEX_ID.index(1)
        end_index = LIST_PATH_MS2.INDEX_ID.index(1000)
        inputListIndex = LIST_PATH_MS2.INDEX_ID[the_index:end_index]
        iID1_list = []
        for i in range(len(inputListIndex)):
            # if inputListEvidence[i] == 0: continue
            iID = inputListIndex[i]
            tmpSeed = self.__captainFillSeed(iID, LIST_PATH_MS2)  # 质荷比，强度
            inputListSeed[i] = tmpSeed
            functionEvidence = CFunctionEvidence()
            # left_mid为MS2的保留时间范围列表中当前MS2左边长度，用来确定当前MS2在保留时间范围内的索引
            tmpEvidence, cycle, index_list_ms2, left_list_ms2, right_list_ms2 = functionEvidence.fillEvidence2(LIST_PATH_MS2, iID, tmpSeed)
            inputListEvidence[i] = tmpEvidence

            iID1 = iID - (iID - 1) % cycle
            if iID1 not in iID1_list:
                iID1_list.append(iID1)
                tmpSeed1 = self.__captainFillSeed(iID1, LIST_PATH_MS1)
                # pro = GetCharge()
                # tmpSeed1.MOZ_FRAG, tmpSeed1.INT_FRAG = pro.get_iso_sin_peak(tmpSeed1)
                inputListSeed1[i] = tmpSeed1
                functionEvidence = CFunctionEvidence()
                tmpEvidence = functionEvidence.fillEvidence1(LIST_PATH_MS1, iID1, tmpSeed1, left_list_ms2, right_list_ms2)
                for j in range(cycle-1):
                    inputListEvidence1[i + j] = tmpEvidence

        return inputListIndex, cycle

    # 将所有的MS2按同一cycle分块
    def cut_list(self, lists, cut_len):
        res_data = []
        if len(lists) > cut_len:
            for i in range(int(len(lists) / cut_len)):
                cut_a = lists[cut_len * i:cut_len * (i + 1)]
                res_data.append(cut_a)
            last_data = lists[int(len(lists) / cut_len) * cut_len:]
            if last_data:
                res_data.append(last_data)
        else:
            res_data.append(lists)
        # yield res_data
        return res_data


    def draw(self):
        listIndex = []

        PathMS2 = CFunctionParseMS2().ms2topkl("D:\\gogogo\\DIA_DeNovo\\DIA\\data\\HeLa-1h_DIA.mzML")
        pathMS2 = CFunctionParseMS2().loadPKL("D:\\gogogo\\DIA_DeNovo\\DIA\\data\\HeLa-1h_DIA.mzML")
        PathMS1 = CFunctionParseMS1().ms1topkl("D:\\gogogo\\DIA_DeNovo\\DIA\\data\\HeLa-1h_DIA.mzML")
        pathMS1 = CFunctionParseMS1().loadPKL("D:\\gogogo\\DIA_DeNovo\\DIA\\data\\HeLa-1h_DIA.mzML")

        listEvidence = [[] * 1] * VALUE_MAX_SCAN  # 在函数中是全局的
        listSeed = [[] * 1] * VALUE_MAX_SCAN
        listEvidence1 = [[] * 1] * VALUE_MAX_SCAN  # 在函数中是全局的
        listSeed1 = [[] * 1] * VALUE_MAX_SCAN

        inputListIndex, cycle = self.__captainGetEvidence(pathMS1, pathMS2, listEvidence1, listSeed1, listEvidence,
                                                          listSeed, listIndex)
        return listEvidence1, listEvidence, cycle


def __fill_ms1_ms2(all_pepmass_charge, num, inputEvidence_list1: list, split_inputEvidence_list2: list):
    time_start = time.time()
    global x_list, matrix_profile_frag, matrix_profile_prec, frag_moz_list, prec_moz_list, x1_list, ms1_index, ms2_rt, ms2_moz, ms2_int, ms1_moz, ms1_int, the_ms1_index, the_ms2_index
    # 1MS1 -- nMS2
    cur = num * 21
    ms1Evidence = inputEvidence_list1[cur]

    base_intensity = np.max(ms1Evidence.MATRIX_PROFILE_PRECURSOR)
    print(base_intensity)
    # if base_intensity > 0.01:
    #     matrix_profile_frag = inputEvidence.MATRIX_PROFILE_FRAG / base_intensity
    # else:
    matrix_profile_prec = ms1Evidence.MATRIX_PROFILE_PRECURSOR

    the_ms1_index = ms1Evidence.THE_MS_INDEX
    ms1_moz = ms1Evidence.PREC_MOZ
    ms1_int = ms1Evidence.PREC_INT

    x1_list = ms1Evidence.LIST_RET_TIME
    print(x1_list)
    prec_moz_list = ms1Evidence.LIST_PREC_MOZ

    # MS1所有母离子的色谱曲线

    for tmp in range(len(split_inputEvidence_list2)):
        ms2Evidence = split_inputEvidence_list2[tmp]
        the_ms2_index = ms2Evidence.THE_MS_INDEX
        left_mid = ms2Evidence.LEFT_MID
        mid_win = ms2Evidence.ACTIVATION_CENTER
        win_size = ms2Evidence.WIN_SIZE
        print("正在拆scan为：" + str(the_ms2_index) + "的MS2")
        if len(ms2Evidence.MATRIX_PROFILE_FRAG) != 0:
            base_intensity = np.max(ms2Evidence.MATRIX_PROFILE_FRAG)
            matrix_profile_frag = ms2Evidence.MATRIX_PROFILE_FRAG

            # ax1 = plt.figure().add_subplot(projection='3d')
            x_list = ms2Evidence.LIST_RET_TIME
            frag_moz_list = ms2Evidence.LIST_FRAG_MOZ
            ms2_rt = ms2Evidence.THE_MS_RT
        else:
            base_intensity = 0

        # 计算相似度，确定电荷状态并写入mgf文件
        CFunctionSimilar().captianGetSimilarCharge(all_pepmass_charge, left_mid, matrix_profile_frag, matrix_profile_prec, frag_moz_list,
                                                       prec_moz_list, x_list, x1_list, ms2_rt, the_ms1_index, ms1_moz,
                                                       ms1_int, the_ms2_index, mid_win, win_size)
    time_end = time.time()
    print("一个cycle的拆谱时间为：" + str(time_end - time_start))


if __name__ == "__main__":
    all_pepmass_charge = {}
    Max_Data = 6000000
    for i in range(Max_Data):
        all_pepmass_charge[i] = []
    # pool = multiprocessing.Pool(processes=4)
    # lock = multiprocessing.Manager().Lock()
    time1_start = time.time()
    test = CFunctionDrawDIAXIC()
    listEvidence1, listEvidence2, cycle = test.draw() # listEvidence1,listEvidence2为所有的MS1 MS2的数据信息，一个listEvidence2对应一个listEvidence1
    time1_end = time.time()
    print("读取所有的MS1和MS2数据信息花费的时间：" + str(time1_end - time1_start))
    inputEvidence_list2 = test.cut_list(listEvidence2, cycle - 1)
    result = []
    time2_start = time.time()
    for i, data in enumerate(inputEvidence_list2):
        __fill_ms1_ms2(all_pepmass_charge, i, listEvidence1, inputEvidence_list2[i])
        # pool.apply_async(__fill_ms1_ms2, (lock, all_pepmass_charge, i, listEvidence1, data),)
    time2_end = time.time()
    print("拆谱时间为：" + str(time2_end - time2_start))
    # pool.close()
    # pool.join()


