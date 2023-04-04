import os
import pickle

from pyteomics import mzml

from MSData import CFileMS1, CFileMS2
from MSOperator import op_INIT_CFILE_MS2, op_INIT_CFILE_MS1
from MSSystem import VALUE_ILLEGAL
from get_MS_data import get_RT_from_rawdata_spectrum, filter_spectrum_ms, get_ION_MOB_from_rawdata_spectrum, \
    get_ISO_WIN_from_rawdata_spectrum, filter_spectrum_ms_dia, get_scan, filter_spectrum_ms_id, \
    get_WIN_from_rawdata_spectrum


def soldierBinarySearch(inputList, start, end, number):

    # check
    if number == inputList[end]:  # 最后一个有问题

        return end

    # find
    while start < end:
        mid = (start + end) // 2
        if number < inputList[mid]:
            end = mid
        elif number > inputList[mid]:
            start = mid + 1
        else:
            return mid

    start = start - 1

    return start

def toolFindIndexFromSortedList0(inputList, start, end, number):

    if number < inputList[start]:
        return VALUE_ILLEGAL  # 没找到必须返回一个非法值

    if number > inputList[end]:
        return VALUE_ILLEGAL  # 没找到必须返回一个非法值

    index = soldierBinarySearch(inputList, start, end, number)

    if inputList[index] != number:

        return VALUE_ILLEGAL  # 没找到必须返回一个非法值

    else:

        return index


def toolFindIndexFromSortedList1(inputList, number):

    return toolFindIndexFromSortedList0(inputList, 0, len(inputList)-1, number)



def toolFindNeighborFromSortedList0(inputList, start, end, number,):

    # if number <= inputList[start]:
    if number <= inputList[start]:
        return start

    if number >= inputList[end]:
        return end

    neighbor = soldierBinarySearch(inputList, start, end, number)  # 通常是落在了左边那个索引上

    disLeft = abs(float(inputList[neighbor]) - float(number))
    disRight = abs(float(inputList[neighbor+1]) - float(number))

    if disLeft < disRight:
        return neighbor
    else:
        return neighbor + 1


def toolFindNeighborFromSortedList1(inputList, number):

    return toolFindNeighborFromSortedList0(inputList, 0, len(inputList)-1, number)  # 长度必须-1，保证是个数字


def toolGetWord(inputString, index, d):
    if inputString[0] != d:
        inputString = d + inputString

    if inputString[-1] != d:
        inputString = inputString + d

    p_d = []

    i = 0
    for c in inputString:

        if c == d:
            p_d.append(i)

        i = i + 1

    result = inputString[p_d[index] + 1:p_d[index + 1]]

    return result

class CFunctionParseMS1:


    def loadPKL(self, pathMS1):

        pathPKL = pathMS1 + '.pkl'
        dataMS1 = CFileMS1()

        f_pkl = open(pathPKL, 'rb')
        dataMS1 = pickle.load(f_pkl)
        f_pkl.close()

        return dataMS1

    def ms1topkl(self, pathMS1):
        path_pkl = pathMS1 + '.pkl'
        dataMS1 = CFileMS1()
        op_INIT_CFILE_MS1(dataMS1)
        with open(pathMS1, 'r')as f:
            data_reader = mzml.MzML(pathMS1)
            ms_level = "ms level"
            # i = 0
            for idx, spectrum in enumerate(data_reader):
                scan = get_scan(spectrum)
                if scan >= 8500 and scan % 22 == 1:
                    RT = get_RT_from_rawdata_spectrum(spectrum)
                    mz_array, intensity_array = filter_spectrum_ms_dia(spectrum)
                    dataMS1.MATRIX_PEAK_MOZ[scan] = mz_array
                    dataMS1.MATRIX_PEAK_INT[scan] = intensity_array
                    dataMS1.INDEX_ID.append(scan)
                    dataMS1.LIST_RET_TIME[scan] = RT
                if scan >= 11000: break
                    # i += 1
                    # if i >= 10000: break
                    # if scan == 4995:
                    #     with open('result.mgf', "a") as f:
                    #         data = zip(mz_array, intensity_array)
                    #         sorted_data = sorted(data, key=lambda x: x[0])
                    #         for f1, f2 in sorted_data:
                    #             print(round(f1,5), round(f2,1), file=f)

        f.close()
        f_pkl = open(path_pkl, 'wb')
        pickle.dump(dataMS1, f_pkl)
        f_pkl.close()
        # path_pkl = pathMS1 + '.pkl'
        # dataMS1 = CFileMS1()
        # op_INIT_CFILE_MS1(dataMS1)
        # i = 0
        # with open(pathMS1, 'r')as f:
        #     data_reader = mzml.MzML(pathMS1)
        #     ms_level = "ms level"
        #     for idx, spectrum in enumerate(data_reader):
        #         id = filter_spectrum_ms_id(spectrum)
        #         if id > 38000 and id < 63000:
        #             if spectrum[ms_level] == 1:
        #                 RT = get_RT_from_rawdata_spectrum(spectrum)
        #                 mz_array, intensity_array, ion_lower, ion_upper = filter_spectrum_ms(spectrum, mz_min, mz_max)
        #                 ion_mob = get_ION_MOB_from_rawdata_spectrum(spectrum)
        #                 dataMS1.MATRIX_PEAK_MOZ[id] = mz_array
        #                 dataMS1.MATRIX_PEAK_INT[id] = intensity_array
        #                 dataMS1.ION_MOBILITY[id] = ion_mob
        #                 dataMS1.INDEX_ID.append(id)
        #                 dataMS1.LIST_RET_TIME[id] = RT / 60
        #         #         i += 1
        #         elif id >= 63000:
        #             break
        # f_pkl = open(path_pkl, 'wb')
        # pickle.dump(dataMS1, f_pkl)
        # f_pkl.close()

class CFunctionParseMS2:


    def loadPKL(self, pathMS2):

        pathPKL = pathMS2 + '.pkl'
        dataMS2 = CFileMS2()

        f_pkl = open(pathPKL, 'rb')
        dataMS2 = pickle.load(f_pkl)

        return dataMS2

    # def ms2topkl(self, pathMS2, mz_min, mz_max):
    #     path_pkl = pathMS2 + '.pkl'
    #     dataMS2 = CFileMS2()
    #     op_INIT_CFILE_MS2(dataMS2)
    #     with open(pathMS2, 'r')as f:
    #         data_reader = mzml.MzML(pathMS2)
    #         ms_level = "ms level"
    #
    #         for idx, spectrum in enumerate(data_reader):
    #             id = filter_spectrum_ms_id(spectrum)
    #             if id > 38000 and id < 63000:
    #                 if spectrum[ms_level] == 2:
    #                     RT = get_RT_from_rawdata_spectrum(spectrum)
    #                     mz_array, intensity_array, ion_lower, ion_upper = filter_spectrum_ms(spectrum, mz_min, mz_max)
    #                     ion_mob = get_ION_MOB_from_rawdata_spectrum(spectrum)
    #                     mid_win = get_ISO_WIN_from_rawdata_spectrum(spectrum)
    #                     id1 = int(id) - int(id) % 65
    #                     dataMS2.MATRIX_PEAK_MOZ[id] = mz_array
    #                     dataMS2.MATRIX_PEAK_INT[id] = intensity_array
    #                     dataMS2.ION_MOBILITY[id] = ion_mob
    #                     dataMS2.INDEX_ID.append(id)
    #                     dataMS2.LIST_ACTIVATION_CENTER[id] = mid_win
    #                     dataMS2.INDEX_RT.append(RT / 60)
    #                     dataMS2.LIST_RET_TIME[id] = RT / 60
    #                     dataMS2.LIST_PRECURSOR_ID.append(id1)
    #                 # i += 1
    #                 elif id >= 63000:
    #                     break
    #     f_pkl = open(path_pkl, 'wb')
    #     pickle.dump(dataMS2, f_pkl)
    #     f_pkl.close()

    # 不含离子淌度的dia数据读取
    def ms2topkl(self, pathMS2):
        path_pkl = pathMS2 + '.pkl'
        dataMS2 = CFileMS2()
        op_INIT_CFILE_MS2(dataMS2)
        with open(pathMS2, 'r')as f:
            data_reader = mzml.MzML(pathMS2)
            ms_level = "ms level"
            # i = 0
            for idx, spectrum in enumerate(data_reader):
                scan = get_scan(spectrum)
                if scan >= 8500 and scan % 22 != 1:
                    RT = get_RT_from_rawdata_spectrum(spectrum)
                    mid_win, win_size = get_WIN_from_rawdata_spectrum(spectrum)
                    dataMS2.LIST_RET_TIME[scan] = RT
                    mz_array, intensity_array = filter_spectrum_ms_dia(spectrum)
                    dataMS2.MATRIX_PEAK_MOZ[scan] = mz_array
                    dataMS2.MATRIX_PEAK_INT[scan] = intensity_array
                    dataMS2.INDEX_ID.append(scan)
                    dataMS2.INDEX_RT.append(RT)
                    dataMS2.WIN_SIZE[scan] = win_size
                    dataMS2.ACTIVATION_CENTER[scan] = mid_win
                if scan >= 11000:break
                    # i += 1
                    # if i >= 10000: break
        f.close()
        f_pkl = open(path_pkl, 'wb')
        pickle.dump(dataMS2, f_pkl)
        f_pkl.close()


