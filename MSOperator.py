from MSData import CFileMS2, CFileMS1
from MSSystem import VALUE_ILLEGAL, VALUE_MAX_SCAN


def op_INIT_CFILE_MS1(dataMS1: CFileMS1):
    dataMS1.INDEX_ID = []
    dataMS1.INDEX_RT = []
    dataMS1.LIST_RET_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    dataMS1.ION_MOBILITY = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    dataMS1.MATRIX_PEAK_INT = [[] * 1] * VALUE_MAX_SCAN
    dataMS1.MATRIX_PEAK_MOZ = [[] * 1] * VALUE_MAX_SCAN


def op_INIT_CFILE_MS2(dataMS2: CFileMS2):
    dataMS2.INDEX_ID= [] # dia数据存的是scan号，从1开始
    dataMS2.INDEX_RT = []
    dataMS2.ION_MOBILITY = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    dataMS2.ACTIVATION_CENTER = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    dataMS2.WIN_SIZE = [VALUE_ILLEGAL] * VALUE_MAX_SCAN

    dataMS2.LIST_RET_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    dataMS2.LIST_PRECURSOR_ID = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    # dataMS2.LIST_ACTIVATION_CENTER = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    # arry = []  # 创建一个空列表
    # for i in range(VALUE_MAX_SCAN):
    #     arry.append([])  # 在空列表中再添加一个空列表
    # dataMS2.MATRIX_PEAK_MOZ = arry
    # arry1 = []  # 创建一个空列表
    # for i in range(VALUE_MAX_SCAN):
    #     arry1.append([])
    # dataMS2.MATRIX_PEAK_INT = arry1
    dataMS2.MATRIX_PEAK_MOZ = [[] * 1] * VALUE_MAX_SCAN
    dataMS2.MATRIX_PEAK_INT = [[] * 1] * VALUE_MAX_SCAN

def opGetStartAndEndForProfile(input_profile, input_seed, input_cutoff, input_n_hole):  # 就这个名字格式特殊

    nHoleLeft = 0
    nHoleRight = 0

    i_middle = input_seed

    if i_middle > 0:
        i_left = i_middle - 1
    else:
        i_left = i_middle

    i_right = i_middle

    result = [i_left, i_right]  # start and end

    int_left = input_profile[i_left]  # int ->当前左边的强度
    int_right = input_profile[i_right]

    int_max = int_left
    int_thr = int_max * input_cutoff

    walkLeft = True
    walkRight = True

    while True:

        if walkLeft or walkRight:
            pass
        else:
            break

        if i_left > 0:

            if walkLeft:
                i_left = i_left - 1

        else:

            walkLeft = False

        if i_right < len(input_profile) - 1:

            if walkRight:
                i_right = i_right + 1

        else:

            walkRight = False

        int_left = input_profile[i_left]
        int_right = input_profile[i_right]

        # max
        if int_max < int_left:
            int_max = int_left

        if int_max < int_right:
            int_max = int_right

        int_thr = int_max * input_cutoff

        # hole
        if int_left < int_thr and walkLeft:
            nHoleLeft = nHoleLeft + 1

        if int_right < int_thr and walkRight:
            nHoleRight = nHoleRight + 1

        if nHoleLeft == input_n_hole:

            walkLeft = False

        if nHoleRight == input_n_hole:

            walkRight = False

    result[0] = i_left
    result[1] = i_right

    return result
