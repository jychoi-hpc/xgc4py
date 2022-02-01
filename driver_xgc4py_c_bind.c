#include <Python.h>
#include "xgc4py_c_bind.h"

int
main() 
{
    PyImport_AppendInittab("xgc4py_c_bind", PyInit_xgc4py_c_bind);
    Py_Initialize();
    PyImport_ImportModule("xgc4py_c_bind");
    call_hello(77);
    xgc4py_init("d3d_coarse_v2", 420);

    double data[6] = {1.,2.,3.,4.,5.,6.};
    long shape[2] = {2,3};
    xgc4py_test_print(data, 3);
    xgc4py_test_print2(data, shape, 2);

    /*
    den, upara, Tperp, Tpara, fn0, fT0 = \
        xgcexp.f0_diag_future(f0_inode1=f0_inode1, ndata=ndata, isp=1, f0_f=f0_f, progress=True)
    print (den.shape, upara.shape, Tperp.shape, Tpara.shape, fn0.shape, fT0.shape)
    (16395, 39, 39) (16395, 39, 39) (16395, 39, 39) (16395, 39, 39) (16395,) (16395,)
     */
    int f0_f_offset = 0;
    int f0_f_ndata = 16395;
    int isp = 1;
    double *f0_f;
    long f0_f_shap[3] = {16395, 39, 39};
    int f0_f_ndim = 3;

    f0_f = malloc(16395*39*39*sizeof(double));
    for (int i=0; i<16395*39*39; i++)
    {
        f0_f[i] = 1.0;
    }

    double *den = malloc(16395*39*39*sizeof(double));

    xgc4py_f0_diag(f0_f_offset, f0_f_ndata, isp, f0_f, f0_f_shap, f0_f_ndim, den);

    Py_Finalize();
    return 0;
}