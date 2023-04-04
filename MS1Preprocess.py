# -*- coding: utf-8 -*-
from MSData import CINI
import numpy as np
import math
# 对MS1预处理：去同位素峰簇转成单同位素峰
# 和确认电荷状态部分合并


class GetCharge:

    def __init__(self):
        self.msms_tol = 20
        self.msms_ppm = 1

        if self.msms_ppm == 1:
            self.msms_fraction = self.msms_tol / 1e6
        else:
            self.msms_fraction = 1

    def get_iso_sin_peak(self, dataMS1Spectrum):
        global charge

        delete_index_list = []
        INI = CINI()
        max_charge = INI.MAX_CHARGE  # 当前考虑1、2、3、4电荷
        isotopic_tol = [INI.MASS_NEUTRON_AVRG / c for c in range(1, max_charge + 1)]
        charge_state = -1
        charge_moz_dic = {}
        the_moz_list = list(dataMS1Spectrum.MOZ_FRAG)
        the_int_list = list(dataMS1Spectrum.INT_FRAG)
        new_peak_int = []
        new_peak_moz = []
        
        while 0 < len(the_moz_list):
            cluster_index_list = []
            # 得到最高峰信号地址
            max_int = -1
            max_int_peak_index = -1
            for i in range(len(the_moz_list)):
                if the_int_list[i] > max_int:
                    max_int = the_int_list[i]
                    max_int_peak_index = i
            # 左右检测是否有同位素峰簇
            # -------------------------- right -----------------------------
            # ##############################################################
            tmp_check_index = max_int_peak_index
            # index = 0
            cluster_index_list.append(tmp_check_index)
            cluster_tail_moz = the_moz_list[tmp_check_index]
            tmp_check_index += 1
            # charge == -1: cluster is not complete

            while tmp_check_index < len(the_moz_list):
                print("finding right cluster charge state+++++++")
                peak_tolerance = the_moz_list[tmp_check_index] - cluster_tail_moz
                if self.msms_ppm == 1:
                    ppm_tol_da = the_moz_list[tmp_check_index] * self.msms_fraction
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
                                cluster_tail_moz = the_moz_list[tmp_check_index]
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
                            nearest_int = the_int_list[cluster_index_list[-1]]
                            if the_int_list[tmp_check_index] > 1.2 * nearest_int:
                                tmp_check_index += 1
                                continue
                            cluster_index_list.append(tmp_check_index)
                            cluster_tail_moz = the_moz_list[tmp_check_index]
                            tmp_check_index += 1
                        # 超过了枚举电荷的范围，break
                        elif peak_tolerance - isotopic_tol[charge_state - 1] > prep_tolerance:
                            break
                        # 又不匹配，又不越界，啥也不是，继续走吧
                        else:
                            tmp_check_index += 1
            # while index < len(dataMS2Spectrum.MOZ_FRAG) over.

            # -------------------------- left ------------------------------
            # ##############################################################
            tmp_check_index = max_int_peak_index - 1
            cluster_left_moz = the_moz_list[max_int_peak_index]
            while tmp_check_index >= 0:
                print("finding left cluster charge state+++++++")
                peak_tolerance = cluster_left_moz - the_moz_list[tmp_check_index]
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
                                nearest_int = the_int_list[cluster_index_list[0]]
                                # if dataMS2Spectrum.INT_FRAG[tmp_check_index] < 0.2 * nearest_int:
                                if the_int_list[tmp_check_index] < 0.3 * nearest_int:
                                    # break
                                    tmp_check_index -= 1
                                    continue

                                charge_state = tol_index + 1
                                # 确定质荷比信息
                                # cluster_index_list.append(tmp_check_index)
                                cluster_index_list = [tmp_check_index] + cluster_index_list
                                cluster_left_moz = the_moz_list[tmp_check_index]
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
                            nearest_int = the_int_list[cluster_index_list[0]]
                            # if dataMS2Spectrum.INT_FRAG[tmp_check_index] < 0.2 * nearest_int:
                            if the_int_list[tmp_check_index] < 0.4 * nearest_int:
                                # break
                                tmp_check_index -= 1
                                continue
                            cluster_index_list = [tmp_check_index] + cluster_index_list
                            cluster_left_moz = the_moz_list[tmp_check_index]
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
                add_moz = the_moz_list.pop(cluster_index_list[0])
                add_int = the_int_list.pop(cluster_index_list[0])
                new_peak_moz.append(add_moz)
                new_peak_int.append(add_int)
            else:
                add_int = 0
                for buf_index in reversed(cluster_index_list[1:]):
                    try:
                        the_moz_list.pop(buf_index)
                        add_int += the_int_list.pop(buf_index)

                    except:
                        pass
                add_moz = the_moz_list.pop(cluster_index_list[0]) * charge_state
                add_moz -= (charge_state - 1) * CINI.MASS_PROTON_MONO
                add_int += the_int_list.pop(cluster_index_list[0])

                new_peak_moz.append(add_moz)
                new_peak_int.append(add_int)

            charge_state = -1
            # cluster 处理结束
        # while over
        # -------------------------------------------------------------------

        index_order = np.argsort(new_peak_moz)
        new_peak_moz = [new_peak_moz[index] for index in index_order]
        new_peak_int = [new_peak_int[index] for index in index_order]

        out_peak_moz = []
        out_peak_int = []

        add_moz = 0
        add_int = 0
        i = 0
        jump_set = set()
        while i < len(new_peak_moz) - 1:
            add_moz = new_peak_moz[i]
            add_int = new_peak_int[i]

            if i in jump_set:
                i += 1
                continue
            for j in range(i + 1, len(new_peak_moz)):
                if self.msms_ppm == 1:
                    prep_tolerance = new_peak_moz[j] * self.msms_fraction
                else:
                    prep_tolerance = self.msms_tol

                if abs(new_peak_moz[i] - new_peak_moz[j]) < prep_tolerance:
                    add_moz = add_moz * add_int + new_peak_moz[j] * new_peak_int[j]
                    add_int += new_peak_int[j]
                    add_moz /= add_int
                    i = j
                    jump_set.add(j)
                elif abs(new_peak_moz[i] - new_peak_moz[j]) >= prep_tolerance:
                    out_peak_moz.append(add_moz)
                    out_peak_int.append(add_int)
                    i = j
                    break

        if add_moz in out_peak_moz:
            pass
        else:
            out_peak_moz.append(add_moz)
            out_peak_int.append(add_int)

        if len(new_peak_moz) == 0:
            pass
        elif jump_set:
            if max(jump_set) == len(new_peak_moz):
                pass
            else:
                out_peak_moz.append(new_peak_moz[-1])
                out_peak_int.append(new_peak_int[-1])
        else:
            out_peak_moz.append(new_peak_moz[-1])
            out_peak_int.append(new_peak_int[-1])


        the_moz_list = out_peak_moz
        the_int_list = out_peak_int
        return the_moz_list, the_int_list

# if __name__ == "__main__":
#     moz = [200.074, 200.10303, 200.13956, 201.08682, 201.09686, 201.12338, 202.08232, 202.10896, 202.11836, 202.12737, 203.06635, 203.10254, 204.0704, 204.1021, 204.11229, 204.13431, 205.09903]
#     it = [2207.1672, 12839.208, 10924.477, 15743.644, 4264.379, 62808.895, 8796.586, 4522.1436, 39202.668, 3480.47, 86626.34, 9371.574, 3174.313, 3814.0977, 2164.832, 40622.36, 16940.693]
#     ms1 = GetCharge()
#     ms1pro = ms1.get_iso_sin_peak(moz, it)
#     print(ms1pro)
