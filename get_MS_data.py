# 读取mzML文件，得到所需要的MS1数据信息
import bisect

from pyteomics import mzml
from MSData import CFileMS1, CFileMS2
import numpy as np
import pandas as pd
import json

ms2_name = "HeLa.ms2"
ms1_name = "HeLa.ms1"


def filter_spectrum_ms_id(spectrum):
    id = spectrum["id"]
    return int(str(id).split(" ")[0][7:])


def filter_spectrum_ms(spectrum, mz_min, mz_max):
    id = spectrum["id"]
    ion_mob_lower = spectrum["ion mobility lower limit"]
    ion_mob_upper = spectrum["ion mobility upper limit"]
    intensity_array = spectrum['intensity array']
    mz_array  = spectrum['m/z array'][intensity_array > 0]
    intensity_array = intensity_array[intensity_array > 0]
    ms_range = (mz_array >= mz_min) & (mz_array < mz_max)
    mz_array  = mz_array[ms_range]
    intensity_array = intensity_array[ms_range]
    return mz_array, intensity_array, ion_mob_lower, ion_mob_upper

# 不含有ion mobility的数据进行测试
def filter_spectrum_ms_dia(spectrum):
    intensity_array = spectrum['intensity array']
    mz_array = spectrum['m/z array'][intensity_array > 0]
    mz = np.array(mz_array, dtype='float32')
    intensity_array = intensity_array[intensity_array > 0]
    intensity = np.array(intensity_array, dtype='float32')
    # ms_range = (mz_array >= mz_min) & (mz_array < mz_max)
    # mz_array  = mz_array[ms_range]
    # intensity_array = intensity_array[ms_range]
    return mz, intensity

def get_scan(spectrum):
    id = int(spectrum["index"])
    return id + 1

# 电荷
def get_charge_state_spectrum(spectrum):
    if "charge state" in spectrum["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]:
        return spectrum["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]['charge state']

# 保留时间：scan start time
def get_RT_from_rawdata_spectrum(spectrum):
    return spectrum["scanList"]["scan"][0]["scan start time"]

# isolation window
def get_WIN_from_rawdata_spectrum(spectrum):
    mid = spectrum["precursorList"]["precursor"][0]["isolationWindow"]["isolation window target m/z"]
    size = spectrum["precursorList"]["precursor"][0]["isolationWindow"]["isolation window lower offset"]
    return mid, size

# 离子淌度
def get_ION_MOB_from_rawdata_spectrum(spectrum):
    return spectrum["scanList"]["scan"][0]["inverse reduced ion mobility"]

# 4D 碎片离子来自于哪个isolation window，即为target m/z加减12.5
def  get_ISO_WIN_from_rawdata_spectrum(spectrum):
    return spectrum["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]['selected ion m/z']

ms_dict = {}
def get_ms2_data(fd):
    with open(ms2_name, "a") as f:
        data_reader = mzml.MzML(fd)
        ms_level = "ms level"
        for idx, spectrum in enumerate(data_reader):
            if spectrum[ms_level] == 2:
                f.write("BEGIN IONS")
                scan = get_scan(spectrum)
                f.write("\n" + "SCAN=" + str(scan))
                RT = get_RT_from_rawdata_spectrum(spectrum)
                PrecursorScan = scan - (scan - 1) % 22
                f.write("\n" + "PRECURSOR SCAN=" + str(PrecursorScan))
                # native_id, mz_array, intensity_array, ion_lower, ion_upper = filter_spectrum_ms(spectrum, mz_min, mz_max)
                mz_array, intensity_array = filter_spectrum_ms_dia(spectrum)
                # ion_mob = get_ION_MOB_from_rawdata_spectrum(spectrum)
                mid_win = get_ISO_WIN_from_rawdata_spectrum(spectrum)
                # id = str(native_id).split(" ")[0][7:]
                # id1 = int(id) - int(id) % 65
                # f.write("\n" + "Native ID:" + id)
                # f.write("\n" + "MS LEVEL=2")
                f.write("\n" + "LIST_ACTIVATION_CENTER=" + str(mid_win))
                # f.write("\n" + "LIST_PRECURSOR_ID=" + str(id1))
                f.write("\n" + "RETENTION TIME=" + str(format(RT * 60, '.6f')) + "\n")
                # f.write("\n" + "ION MOBILITY=" + str(ion_mob) + "\n")
                data = zip(mz_array, intensity_array)
                sorted_data = sorted(data, key=lambda x: x[0])
                for f1,f2 in sorted_data:
                    print(str(format(f1, '.6f')), str(format(f2, '.4f')), file=f)
                f.write("END IONS" + "\n")
    f.close()

def get_ms1_data(fd):
    with open(ms1_name, "a") as f:
        data_reader = mzml.MzML(fd)
        ms_level = "ms level"
        for idx, spectrum in enumerate(data_reader):
            if spectrum[ms_level] == 1:
                f.write("BEGIN IONS")
                scan = get_scan(spectrum)
                f.write("\n" + "SCAN=" + str(scan))
                RT = get_RT_from_rawdata_spectrum(spectrum)
                mz_array, intensity_array = filter_spectrum_ms_dia(spectrum)
                # native_id, mz_array, intensity_array, ion_lower, ion_upper = filter_spectrum_ms(spectrum, mz_min, mz_max)
                # ion_mob = get_ION_MOB_from_rawdata_spectrum(spectrum)
                # id = str(native_id).split(" ")[0][7:]
                # f.write("\n" + "Native ID:" + id)
                # f.write("\n" + "MS LEVEL=2")
                f.write("\n" + "RETENTION TIME=" + str(format(RT * 60, '.6f')) + "\n")
                # f.write("\n" + "ION MOBILITY=" + str(ion_mob) + "\n")
                data = zip(mz_array, intensity_array)
                sorted_data = sorted(data, key=lambda x: x[0])
                for f1,f2 in sorted_data:
                    print(str(format(f1, '.6f')), str(format(f2, '.4f')), file=f)
                f.write("END IONS" + "\n")
    f.close()

path = "D:\\gogogo\\DIA_DeNovo\\DIA\\data\\HeLa-1h_DIA.mzML"

if __name__ == "__main__":
    get_ms1_data(path)
    get_ms2_data(path)
