# cython: language_level=3
# distutils: language = c++

import numpy as np
#import adios2 as ad2
from tqdm import tqdm

import os
import subprocess
import logging

from math import sqrt, floor, exp
#import torch
import sys

from libc.string cimport memcpy

from xgc4py import XGC

cimport numpy as np
np.import_array()

cdef xgcexp

def log(*args, logtype='debug', sep=' '):
    getattr(logging, logtype)(sep.join(map(str, args)))

cdef public void xgc4py_init(char* expdir, int timestep):
    global xgcexp
    logging.basicConfig (
        level = logging.DEBUG,
        format = '[XGC4PY-C:%(levelname)s] %(message)s')

    log (sys.path)
    log ("XGC", XGC)
    log ('xgc4py_init', expdir, timestep)
    xgcexp = XGC("d3d_coarse_v2", step=timestep)
    log (xgcexp)

cdef np.ndarray double_arr_to_numpy(double* data, long* shape, int ndim):
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

    return x

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

cdef public void xgc4py_f0_diag(int f0_inode1, int ndata, int isp, 
            double* f0_f, long* shape, int ndim,
            double* den, double* u_para, double* T_perp, double* T_para):
    """
    This is a wrapper cython function to call f0_diag from C/C++.
    Parameters:
        f0_inode1 : offset
        ndata     : number of mesh nodes to calculate
        isp       : electron(=0) or ion(=1)
        f0_f      : pointer of f0_f double array in C/C++ side
        shape     : shape of f0_f data
        ndim      : number of dimension of f0_f (asume to be 3)
        den       : (output) pointer of density array in C/C++ side. 
                    Assume the memory is already allocated. 
                    It should be same size of f0_f
        u_para    : (output) pointer of u_para array in C/C++ side. Same assuption with den
        T_perp    : (output) pointer of T_perp array in C/C++ side. Same assuption with den
        T_para    : (output) pointer of T_para array in C/C++ side. Same assuption with den
        

    """
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
    memcpy(buff, f0_f, x.nbytes)

    (den_, u_para_, T_perp_, T_para_, n0_, T0_) = xgcexp.f0_diag(f0_inode1, ndata, isp, x)
    log ('f0_diag output:', den_.shape, u_para_.shape, T_perp_.shape, T_para_.shape, n0_.shape, T0_.shape)

    cdef np.ndarray[np.double_t, ndim=3, mode='c'] out_den_buff
    cdef np.ndarray[np.double_t, ndim=3, mode='c'] out_u_para_buff
    cdef np.ndarray[np.double_t, ndim=3, mode='c'] out_T_perp_buff
    cdef np.ndarray[np.double_t, ndim=3, mode='c'] out_T_para_buff

    out_den_buff = np.ascontiguousarray(den_, dtype=np.double)
    memcpy(den, <double*> out_den_buff.data, den_.nbytes)

    out_u_para_buff = np.ascontiguousarray(u_para_, dtype=np.double)
    memcpy(u_para, <double*> out_u_para_buff.data, u_para_.nbytes)

    out_T_perp_buff = np.ascontiguousarray(T_perp_, dtype=np.double)
    memcpy(T_perp, <double*> out_T_perp_buff.data, T_perp_.nbytes)

    out_T_para_buff = np.ascontiguousarray(T_para_, dtype=np.double)
    memcpy(T_para, <double*> out_T_para_buff.data, T_para_.nbytes)

cdef public void xgc4py_f0_param(int f0_inode1, int ndata, int isp, 
            double* f0_f, long* shape, int ndim,
            double* vol, double* vth, double* vth2, double* vp, double* mu_vol,
            double* ptl_mass, double* sml_e_charge,
            double* f0_grid_vol, double* mu_vp_vol):
    """
    This is a wrapper cython function to call f0_diag from C/C++.
    Parameters:
        f0_inode1 : offset
        ndata     : number of mesh nodes to calculate
        isp       : electron(=0) or ion(=1)
        f0_f      : pointer of f0_f double array in C/C++ side
        shape     : shape of f0_f data
        ndim      : number of dimension of f0_f (asume to be 3)
        vol       : (output) pointer of density array in C/C++ side. 
                    Assume the memory is already allocated. 
                    It should be same size of f0_f
        vth, vth2 : (output) shape :(ndata,)
        vp, mu_vol : (output) shape: (39,)
        ptl_mass, sml_e_charge : (output) double
        f0_grid_vol : (output) shape: (ndata,)
        mu_vp_vol : (output) shape: (39,39)
    """
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
    memcpy(buff, f0_f, x.nbytes)

    (vol_, vth_, vp_, mu_vol_, vth2_, ptl_mass_, sml_e_charge_, f0_grid_vol_, mu_vp_vol_) = xgcexp.f0_param(f0_inode1, ndata, isp, x)
    log ('f0_param output:', vol_.shape, vth_.shape, vp_.shape, mu_vol_.shape, vth2_.shape)
    log ('f0_param output:', ptl_mass_, sml_e_charge_, f0_grid_vol_.shape, mu_vp_vol_.shape)

    cdef np.ndarray[np.double_t, ndim=3, mode='c'] out_vol_buff
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] out_vth_buff
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] out_vth2_buff
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] out_vp_buff
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] out_mu_vol_buff
    cdef np.ndarray[np.double_t, ndim=0, mode='c'] out_ptl_mass_buff
    cdef np.ndarray[np.double_t, ndim=0, mode='c'] out_sml_e_charge_buff
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] out_f0_grid_vol
    cdef np.ndarray[np.double_t, ndim=2, mode='c'] out_mu_vp_vol

    out_vol_buff = np.ascontiguousarray(vol_, dtype=np.double)
    memcpy(vol, <double*> out_vol_buff.data, out_vol_buff.nbytes)

    out_vth_buff = np.ascontiguousarray(vth_, dtype=np.double)
    memcpy(vth, <double*> out_vth_buff.data, out_vth_buff.nbytes)

    out_vth2_buff = np.ascontiguousarray(vth2_, dtype=np.double)
    memcpy(vth2, <double*> out_vth2_buff.data, out_vth2_buff.nbytes)

    out_vp_buff = np.ascontiguousarray(vp_, dtype=np.double)
    memcpy(vp, <double*> out_vp_buff.data, out_vp_buff.nbytes)

    out_mu_vol_buff = np.ascontiguousarray(mu_vol_, dtype=np.double)
    memcpy(mu_vol, <double*> out_mu_vol_buff.data, out_mu_vol_buff.nbytes)

    out_ptl_mass_buff = np.array(ptl_mass_, dtype=np.double)
    memcpy(ptl_mass, <double*> out_ptl_mass_buff.data, out_ptl_mass_buff.nbytes)

    out_sml_e_charge_buff = np.array(sml_e_charge_, dtype=np.double)
    memcpy(sml_e_charge, <double*> out_sml_e_charge_buff.data, out_sml_e_charge_buff.nbytes)

    out_f0_grid_vol = np.ascontiguousarray(f0_grid_vol_, dtype=np.double)
    memcpy(f0_grid_vol, <double*> out_f0_grid_vol.data, out_f0_grid_vol.nbytes)

    out_mu_vp_vol = np.ascontiguousarray(mu_vp_vol_, dtype=np.double)
    memcpy(mu_vp_vol, <double*> out_mu_vp_vol.data, out_mu_vp_vol.nbytes)
