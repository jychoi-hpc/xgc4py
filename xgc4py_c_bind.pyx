# cython: language_level=3
# distutils: language = c++

import numpy as np
import adios2 as ad2
from tqdm import tqdm

import os
import subprocess
import logging

from math import sqrt, floor, exp
import torch

from libc.string cimport memcpy

import sys
sys.path.insert(0, '')

from XGC import XGC, hello

cimport numpy as np
np.import_array()

cdef public void call_hello(int x):
    hello(x)

cdef xgcexp

def log(*args, logtype='debug', sep=' '):
    getattr(logging, logtype)(sep.join(map(str, args)))

cdef public void xgc4py_init(char* expdir, int timestep):
    global xgcexp
    logging.basicConfig (
        level = logging.DEBUG,
        format = '[%(levelname)s] %(message)s')

    log ('xgc4py_init', expdir, timestep)
    xgcexp = XGC("d3d_coarse_v2", step=timestep)
    log (xgcexp)

"""
Test with numpy
"""
cdef public void xgc4py_test_print(double* data, int len):
    global xgcexp
    x = np.zeros(len, dtype=np.double)
    cdef np.ndarray[np.double_t, ndim=1, mode = 'c'] np_buff = np.ascontiguousarray(x, dtype=np.double)
    cdef double* buff = <double*> np_buff.data
    memcpy(buff, data, sizeof(double)*len)
    xgcexp.test_print(x)
