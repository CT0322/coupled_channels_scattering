# -*- coding: utf-8 -*-
# rebuild_corr(dir_path, _output_path, file_suffix)
#############################################################################################
# Author: Ting Chen
# Date: 　2017/2/14
# Function：combine different files corresponding to different configurations into one file
#############################################################################################

import sys
import os
import numpy as np


def rebuild_corr(dir_path, _output_path, file_suffix):

    files = os.listdir(dir_path)
    files = filter(lambda x: file_suffix in x, files)
    matrix = []
    head = None

    for _file in files:
        _file_path = os.path.join(dir_path, _file)
        if head is None:
            f0 = open(_file_path, 'r')
            head = f0.readlines(1)
            head = head[0]
            head = head.split()
            head = map(lambda x: int(x), head)
            head[0] = len(files)
            head = map(lambda x: str(x), head)
            head = ' '.join(head)
          #  print head
          #  head = head[1:]
        mat = np.loadtxt(_file_path, skiprows=1)
        matrix.append(mat)
    matrix = np.concatenate(matrix)
    # matrix = np.concatenate([matrix, np.nan*np.ones([len(matrix), len(head)-matrix.shape[1]])], axis=1)
    # matrix = np.concatenate([np.array(head).reshape((1, len(head))), matrix])
    output_path = os.path.join(_output_path, file_suffix + '.dat')
    np.savetxt(output_path, matrix, fmt=['%d', '%.16e', '%.16e'], header=head,comments='')
    #print matrix



