import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn2_circles
import csv
import numpy as np
import pandas as pd

# 1.将整个数组对半分开，直到只剩一个元素
 # 2.从一个元素开始往回进行对比合并，将合并的内容暂且放在一个空列表中
 # 定义合并merge函数
def merge(L,R):
    # 将两个排过序的两个列表进行合并并排序
    # 分别用于限定L、R的数据减少情况
    i, j = 0,0
    # 用于存放L与R的合并内容
    res = []
    while i < len(L) and j < len(R):
        # 如果左边的数大于右边的数，就将左边的数先添加到res中，再继续比较（合并的R、L已经排过序）
        # 如果如果右边的数大于左边的数，就将右边的数先添加到res中，再继续比较（合并的R、L已经排过序）
        if L[i] <= R[j]:
            res.append(L[i])
            i += 1
        else:
            res.append(R[j])
            j += 1
    # 因为当 i == len(L) 或者 j == len(R)时，跳出while循环，且每次循环只处理一个列表里的内容，所以其中有一个列表内容会先全部加入res中，另一个列表还剩内容未加进res中。
    # 对未处理完的列表，直接加入res中
    res += R[j:] if i == len(L) else L[i:]
    return res

 # 定义排序merge_sort函数
def merge_sort(List):
    length = len(List)
    if length <= 1:
        return List
    else:
        mid = length//2
        left = merge_sort(List[:mid])
        right = merge_sort(List[mid:])
        return merge(left,right)






