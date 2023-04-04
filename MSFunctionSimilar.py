# 计算母离子与碎片离子的相似度，做拆分（如何调用DIA-Umpire）
import heapq
import math
import random
import time

import seaborn as sns
import numpy as np
import torch
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from MSData import CINI
from MSSystem import PLOT_COLOR, DRAW_COLOR
import xlwt


class CFunctionGetCharge:

    def __init__(self):
        self.msms_tol = 20
        self.msms_ppm = 1

        if self.msms_ppm == 1:
            self.msms_fraction = self.msms_tol / 1e6
        else:
            self.msms_fraction = 1

    def get_iso_peak(self, pre_moz, ms1_moz, ms1_int):
        global charge

        INI = CINI()
        the_ms1_moz = list(ms1_moz)  # 当前ms1谱图的质荷比
        the_ms1_int = list(ms1_int)  # 当前ms1谱图的谱峰强度
        pre_idx = the_ms1_moz.index(pre_moz)
        max_charge = INI.MAX_CHARGE  # 当前考虑1、2、3、4电荷
        isotopic_tol = [INI.MASS_NEUTRON_AVRG / c for c in range(1, max_charge + 1)]
        charge_state = -1
        charge_moz_dic = {}
        while 0 < len(the_ms1_moz):
            cluster_index_list = []
            # 得到最高峰信号地址
            max_int = -1
            max_int_peak_index = -1
            for i in range(len(the_ms1_moz)):
                if the_ms1_int[i] > max_int:
                    max_int = the_ms1_int[i]
                    max_int_peak_index = i
            # 左右检测是否有同位素峰簇
            # -------------------------- right -----------------------------
            # ##############################################################
            tmp_check_index = max_int_peak_index
            # index = 0
            cluster_index_list.append(tmp_check_index)
            cluster_tail_moz = the_ms1_moz[tmp_check_index]
            tmp_check_index += 1
            # charge == -1: cluster is not complete

            while tmp_check_index < len(the_ms1_moz):
                print("finding right cluster charge state+++++++")
                peak_tolerance = the_ms1_moz[tmp_check_index] - cluster_tail_moz
                if self.msms_ppm == 1:
                    ppm_tol_da = the_ms1_moz[tmp_check_index] * self.msms_fraction
                else:
                    ppm_tol_da = self.msms_tol
                prep_tolerance = ppm_tol_da
                if peak_tolerance < isotopic_tol[-1] - prep_tolerance:
                    # 小于最小,跳过去
                    # 要继续找，不要break
                    tmp_check_index += 1

                elif peak_tolerance > isotopic_tol[0] + prep_tolerance:
                    # 大于最大
                    # 不要继续找，要break
                    break

                else:

                    # 先确定电荷状态，再继续cluster构建
                    if charge_state == -1:
                        for tol_index in range(len(isotopic_tol)):
                            if math.isclose(peak_tolerance, isotopic_tol[tol_index], abs_tol=prep_tolerance):
                                charge_state = tol_index + 1

                                # 确定质荷比信息
                                cluster_index_list.append(tmp_check_index)
                                cluster_tail_moz = the_ms1_moz[tmp_check_index]
                                break
                        if charge_state == -1:
                            # 继续向后寻找（MUST）
                            tmp_check_index += 1
                        else:

                            # 连续向后构造
                            tmp_check_index += 1


                    # 已经确定电荷状态
                    else:
                        # 仍然是我的同位素峰，进入峰簇下标列表
                        if math.isclose(peak_tolerance, isotopic_tol[charge_state - 1], abs_tol=prep_tolerance):

                            # 尽可能避免误匹配的发生,混合谱峰时，拒绝构造高低错落类型谱峰
                            nearest_int = the_ms1_int[cluster_index_list[-1]]
                            if the_ms1_int[tmp_check_index] > 1.2 * nearest_int:
                                tmp_check_index += 1
                                continue
                            cluster_index_list.append(tmp_check_index)
                            cluster_tail_moz = the_ms1_moz[tmp_check_index]
                            tmp_check_index += 1
                        # 超过了枚举电荷的范围，break
                        elif peak_tolerance - isotopic_tol[charge_state - 1] > prep_tolerance:
                            break
                        # 又不匹配，又不越界，啥也不是，继续走吧
                        else:
                            tmp_check_index += 1
            # while index < len(dataMS2Spectrum.LIST_PEAK_MOZ) over.

            # -------------------------- left ------------------------------
            # ##############################################################
            tmp_check_index = max_int_peak_index - 1
            cluster_left_moz = the_ms1_moz[max_int_peak_index]
            while tmp_check_index >= 0 and the_ms1_moz:
                print("finding left cluster charge state+++++++")
                peak_tolerance = cluster_left_moz - the_ms1_moz[tmp_check_index]
                if self.msms_ppm == 1:
                    ppm_tol_da = cluster_left_moz * self.msms_fraction
                else:
                    ppm_tol_da = self.msms_tol
                prep_tolerance = ppm_tol_da
                if peak_tolerance < isotopic_tol[-1] - prep_tolerance:
                    # 小于最小,跳过去
                    # 要继续找，不要break
                    tmp_check_index -= 1

                elif peak_tolerance > isotopic_tol[0] + prep_tolerance:
                    # 大于最大
                    # 不要继续找，要break
                    break

                else:

                    # 先确定电荷状态，再继续cluster构建
                    if charge_state == -1:
                        for tol_index in range(len(isotopic_tol)):
                            if math.isclose(peak_tolerance, isotopic_tol[tol_index], abs_tol=prep_tolerance):

                                # 尽可能避免误匹配的发生,混合谱峰时，拒绝构造高低错落类型谱峰
                                # nearest_int: cluster list中最左端的谱峰信号的强度值
                                nearest_int = the_ms1_int[cluster_index_list[0]]
                                # if dataMS2Spectrum.LIST_PEAK_INT[tmp_check_index] < 0.2 * nearest_int:
                                if the_ms1_int[tmp_check_index] < 0.3 * nearest_int:
                                    # break
                                    tmp_check_index -= 1
                                    continue

                                charge_state = tol_index + 1
                                # 确定质荷比信息
                                # cluster_index_list.append(tmp_check_index)
                                cluster_index_list = [tmp_check_index] + cluster_index_list
                                cluster_left_moz = the_ms1_moz[tmp_check_index]
                                break
                        if charge_state == -1:
                            # 继续向前寻找（MUST）
                            tmp_check_index -= 1
                        else:

                            # 连续向后构造
                            tmp_check_index -= 1


                    # 已经确定电荷状态
                    else:
                        # 仍然是我的同位素峰，进入峰簇下标列表
                        if math.isclose(peak_tolerance, isotopic_tol[charge_state - 1], abs_tol=prep_tolerance):
                            # 尽可能避免误匹配的发生,混合谱峰时，拒绝构造高低错落类型谱峰
                            nearest_int = the_ms1_int[cluster_index_list[0]]
                            # if dataMS2Spectrum.LIST_PEAK_INT[tmp_check_index] < 0.2 * nearest_int:
                            if the_ms1_int[tmp_check_index] < 0.4 * nearest_int:
                                # break
                                tmp_check_index -= 1
                                continue
                            cluster_index_list = [tmp_check_index] + cluster_index_list
                            cluster_left_moz = the_ms1_moz[tmp_check_index]
                            tmp_check_index -= 1
                        # 超过了枚举电荷的范围，break
                        elif peak_tolerance - isotopic_tol[charge_state - 1] > prep_tolerance:
                            break
                        # 又不匹配，又不越界，啥也不是，继续走吧
                        else:
                            tmp_check_index -= 1

            # --------------------------------------------------------------
            # cluster 构造结束，开始收工
            print("charge state: ", charge_state)
            if charge_state == -1:
                # 删除
                # 没找到电荷状态默认为2+
                the_ms1_moz.pop(cluster_index_list[0])
                the_ms1_int.pop(cluster_index_list[0])
            else:
                if charge_state not in charge_moz_dic.keys():
                    charge_moz_dic[charge_state] = cluster_index_list
                else:
                    for i in cluster_index_list:
                        charge_moz_dic[i] = charge_state
                    sort_list = sorted(cluster_index_list, reverse=True)
                    for idx in sort_list:
                        the_ms1_moz.pop(idx)
                        the_ms1_int.pop(idx)

            charge_state = -1
            # cluster 处理结束
        # while over
        print(charge_moz_dic)
        # for key, value in charge_moz_dic.items():
        #     if pre_idx in value:
        #         charge = key
        #         return charge
        # charge = 2
        # return charge
        charge = charge_moz_dic.get(pre_idx)
        if charge is None:
            charge = 2
        return charge
        # -------------------------------------------------------------------


def randomcolor():
    colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    color = ""
    for i in range(6):
        color += colorArr[random.randint(0,14)]
    return "#"+color

def mean_processed(a):
    b = np.zeros_like(a)
    num = 0
    for i, j in enumerate(a):
        if j.all() != 0:
            num += 1

    for i, j in enumerate(a):
        if j != 0:
            b[i] = np.sum(a) / num
    #     print(np.sum(a) / num)
    #     print(b)
    return b


def Pearson_similar(item1, item2):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    result = [[] * 1] * len(item1)
    for j, pre in enumerate(item1):
        tmp_list = []
        for i, xic in enumerate(item2):
            x = torch.tensor(np.array(pre), dtype=torch.float64, device=device)
            y = torch.tensor(np.array(xic), dtype=torch.float64, device=device)
            cos = torch.nn.CosineSimilarity(dim = 0)
            res = cos(x, y)
            similar = res
            similar.cpu()
            tmp_list.append(format(similar, '.2f'))
            result[j] = tmp_list
            # pre = pre.astype(np.float)
            # xic = xic.astype(np.float)
            # x = torch.tensor(np.array([pre, xic]), dtype=torch.float64, device=device)
            # res = torch.corrcoef(x)
            # similar = res[0][1]
            # similar.float()
            # similar.cpu()
            # tmp_list.append(format(similar, '.2f'))
            # result[j] = tmp_list
    return result


def sum_pearson(result_list):
    j = 0
    for i, data in enumerate(result_list):
        if float(data) > 0.6:
            j += 1
    return j

def corr_index(result_list):
    index = []
    num = 0
    for i, data in enumerate(result_list):
        if float(data) > 0.6:
            index.append(i)
            num += 1
    return index, num


class CFunctionSimilar:

    def print_mgf(self, charge, flag, the_ms2_index, the_ms2_moz, the_ms2_int, pre_moz, pre_int, the_ms2_rt):
        # 最初1电荷2电荷3电荷都为0，如果有一个1+电荷的生成mgf了，则flag1+1
        with open('10000-10500.mgf', 'a') as f:
            f.write('BEGIN IONS\n')
            # f.write('TITLE=HeLa.' + str(the_ms2_index) + '.' + str(the_ms2_index) + '.' + str(charge) + '.' + str(flag) + " " + "File:\"HeLa-1h_DIA.raw\"" + "," + " NativeID:\"controllerType=0 controllerNumber=1 " + "scan=" + str(the_ms2_index) + "\"" + '\n')
            f.write('TITLE=HeLa.' + str(the_ms2_index) + '.' + str(the_ms2_index) + '.' + str(charge) + '.' + str(flag) + '.dta' + '\n')
            f.write('CHARGE=' + str(charge) + '+' + '\n')
            f.write('RTINSECONDS=' + str(format(the_ms2_rt * 60, '.7f')) + '\n')
            # f.write('PEPMASS=' + str(format(pre_moz, '.6f')) + " " + str(format(pre_int, '.4f')) + '\n')
            f.write('PEPMASS=' + str(format(pre_moz, '.6f')) + '\n')
            data = zip(the_ms2_moz, the_ms2_int)
            for f1, f2 in data:
                print(str(format(f1, '.5f')), str(format(f2,'.1f')), file=f)
            f.write('END IONS' + '\n')

    def captianGetSimilarCharge(self, all_pepmass_charge, left_mid, frag_profile, all_prec_profile, frag_moz_list, all_prec_moz_list, frag_ret, prec_ret, the_ms2_rt,
                                  the_ms1_index, the_ms1_moz, the_ms1_int, the_ms2_index, mid_win, win_size):
        global charge, ms2_moz, ms2_int
        index = []
        # 获取当前MS2对应的isolation window，得到母离子的质荷比范围
        for i, pre_mz in enumerate(all_prec_moz_list):
            if pre_mz[left_mid] >= mid_win - win_size and pre_mz[left_mid] <= mid_win + win_size:
                index.append(i)
        prec_profile = np.zeros(shape=[len(index), len(prec_ret)])
        prec_moz_list = np.zeros(shape=[len(index), len(prec_ret)])
        for j, cur in enumerate(index):
            prec_moz_list[j] = all_prec_moz_list[cur]
            prec_profile[j] = all_prec_profile[cur]
        tmp_list = [[] * 1] * len(prec_moz_list)
        start = time.time()
        # result_list = Pearson_similar(prec_profile, frag_profile)
        result_list = Pearson_similar(all_prec_profile, frag_profile)
        end = time.time()
        print("比较时间花费：" + str(end - start))
        # 取相关系数大于0.6的碎片离子
        num_list = []  # 记录每一根母离子与碎片离子计算相关性，碎片离子相关系数大于0.6的谱峰个数
        for i, cor in enumerate(result_list):
            tmp, num = corr_index(result_list[i])  # tmp为大于0.6的碎片离子索引， num为result[i]这一行中大于0.6的个数
            tmp_list[i] = tmp
            num_list.append(num)
        del result_list
        # ---------------------------------------------------------------------------------
        sort_num = sorted(enumerate(num_list), key=lambda num_list: num_list[1],
                          reverse=False)  # 按相关性从小到大排序 sort_num[0]为索引值，[1]为相关性
        sort_index = [num_list[0] for num_list in sort_num]  # 相关性从小到大排序对应的母离子下标
        pep_mass_list = {}  # key: pepmass value: 相关性
        for i in range(len(sort_index)):
            pep_mass = prec_moz_list[sort_index[i]][left_mid]
            pep_mass_list[pep_mass] = [pep_mass, sort_num[i][1]]  # 存放的是每一根PEPMASS以及它的相关性[825.4846801757812, [825.4846801757812, 0]]
        # 先对所有可能的PEPMASS去冗余，首先谱图内部去冗余，然后谱图之间去冗余
        # 谱图内部去冗余
        new_pep_mass_list = {}
        for k, v in pep_mass_list.items():
            k = '{:.3f}'.format((float(k)))
            new_pep_mass_list[k] = v  # v[0]存放的即为去冗余之后的pepmass
        sort_pep_mass = []
        sort_pep_mass_sim = []
        for data in new_pep_mass_list.values():
            sort_pep_mass.append(data[0])
            sort_pep_mass_sim.append(data[1])
        sort_pep_mass = list(reversed(sort_pep_mass))  # 翻转 相关性从大到小的PEPMASS
        # 接下来是谱图之间的去冗余，首先需要获取每一个pepmass的电荷状态
        pepmass_charge = {}  # 2+
        # pepmass_charge2 = {}  # 3+
        for i, mz in enumerate(sort_pep_mass):
            charge = CFunctionGetCharge().get_iso_peak(mz, the_ms1_moz, the_ms1_int)
            pepmass_charge[sort_pep_mass[i]] = charge
            # pepmass_charge2[sort_pep_mass[i]] = charge2
        result_pep_mass = []

        # 去all集合里查找是否重复
        for k, v in pepmass_charge.items():
            label = False  # 标志 True表示重复
            key = int(float(k) * 100)
            if (len(result_pep_mass) == 30): break  # 不重复的topk
            for i in range(key - 2, key + 3):
                if v in all_pepmass_charge[i]:
                    label = True
                    break
            if label == False:
                all_pepmass_charge[key].append(v)
                result_pep_mass.append(k)
        # f1 = open("pepmass.txt", "a+")
        # f1.writelines(str(result_pep_mass))

        tmp_col = list(prec_moz_list[:, left_mid])
        sort_index = []

        for i, data in enumerate(result_pep_mass):
            sort_index.append(tmp_col.index(data))

        # -----------------------------------------------------------------

        flag = 0

        for i, index in enumerate(sort_index):
            tmp = tmp_list[index]
            ms2_moz = []
            ms2_int = []
            pre_moz = prec_moz_list[index][left_mid]
            pre_int = prec_profile[index][left_mid]
            # 获取这跟母离子对应的相关性>0.6的碎片离子
            for j, data in enumerate(tmp):
                ms2_int.append(frag_profile[data][left_mid])
                ms2_moz.append(frag_moz_list[data][left_mid])
            ms2_data = (zip(ms2_moz, ms2_int))
            sorted_data = sorted(ms2_data, key=lambda x: x[0])
            the_ms2_moz = []
            the_ms2_int = []
            # the_ms2 存放的是当前MS2的质荷比和谱峰强度
            # 去除特殊情况
            for f1, f2 in sorted_data:
                if f1 == 0.00000: continue
                the_ms2_moz.append(f1)
                the_ms2_int.append(f2)

            if the_ms2_moz == [] and the_ms2_int == []: return
            charge = pepmass_charge[pre_moz]  # 获取当前母离子的电荷状态
            CFunctionSimilar().print_mgf(charge, flag, the_ms2_index, the_ms2_moz, the_ms2_int, pre_moz, pre_int, the_ms2_rt)
            flag += 1

