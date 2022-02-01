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

cdef public void xgc4py_test_print2(double* data, long* shape, int ndim):
    global xgcexp

    ## shape setup
    z = np.zeros(ndim, dtype=np.int)
    cdef np.ndarray[np.int_t, ndim=1, mode='c'] z_buff = np.ascontiguousarray(z, dtype=np.int)
    cdef long* z_buff_data = <long*> z_buff.data
    memcpy(z_buff_data, shape, z.nbytes)

    ## f0 setup 
    ## This is a workaround to set ndim
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] x_buff1
    cdef np.ndarray[np.double_t, ndim=2, mode='c'] x_buff2
    cdef np.ndarray[np.double_t, ndim=3, mode='c'] x_buff3
    cdef np.ndarray[np.double_t, ndim=4, mode='c'] x_buff4
    cdef np.ndarray[np.double_t, ndim=5, mode='c'] x_buff5
    cdef double* buff

    x = np.zeros(z, dtype=np.double)
    if ndim == 1:
        x_buff1 = np.ascontiguousarray(x, dtype=np.double)
        buff = <double*> x_buff1.data
    if ndim == 2:
        x_buff2 = np.ascontiguousarray(x, dtype=np.double)
        buff = <double*> x_buff2.data
    elif ndim == 3:
        x_buff3 = np.ascontiguousarray(x, dtype=np.double)
        buff = <double*> x_buff3.data
    elif ndim == 4:
        x_buff4 = np.ascontiguousarray(x, dtype=np.double)
        buff = <double*> x_buff4.data
    elif ndim == 5:
        x_buff5 = np.ascontiguousarray(x, dtype=np.double)
        buff = <double*> x_buff5.data
    memcpy(buff, data, x.nbytes)

    xgcexp.test_print(z)
    xgcexp.test_print(x)
