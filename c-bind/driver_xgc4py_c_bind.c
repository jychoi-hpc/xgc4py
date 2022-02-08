#include <Python.h>
#include "xgc4py_c_bind.h"

#include <adios2.h>

int
main() 
{
    /* Setup */
    PyImport_AppendInittab("xgc4py_c_bind", PyInit_xgc4py_c_bind);
    Py_Initialize();
    PyImport_ImportModule("xgc4py_c_bind");

    /* Initialization */
    xgc4py_init("d3d_coarse_v2", 420);

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

    f0_f = (double*) malloc(16395*39*39*sizeof(double));
    for (int i=0; i<16395*39*39; i++)
    {
        f0_f[i] = 1.0;
    }

    double *vol = (double*) malloc(16395*39*39*sizeof(double));
    double *vth = (double*) malloc(16395*sizeof(double));
    double *vth2 = (double*) malloc(16395*sizeof(double));
    double *vp = (double*) malloc(39*sizeof(double));
    double *mu_vol = (double*) malloc(39*sizeof(double));
    double ptl_mass, sml_e_charge;
    double *f0_grid_vol = (double*) malloc(16395*sizeof(double));
    double *mu_vp_vol = (double*) malloc(39*39*sizeof(double));

    /* Call f0_param */
    xgc4py_f0_param(f0_f_offset, f0_f_ndata, isp, 
                f0_f, f0_f_shap, f0_f_ndim,  /* f0_f data (input) */
                vol, vth, vth2, vp, mu_vol, /* (output) */
                &ptl_mass, &sml_e_charge, /* (output) */
                f0_grid_vol, mu_vp_vol); /* (output) */
    printf ("xgc4py_f0_param: ptl_mass= %g\n", ptl_mass);
    printf ("xgc4py_f0_param: sml_e_charge= %g\n", sml_e_charge);

    double *den = (double*) malloc(16395*39*39*sizeof(double));
    double *u_para = (double*) malloc(16395*39*39*sizeof(double));
    double *T_perp = (double*) malloc(16395*39*39*sizeof(double));
    double *T_para = (double*) malloc(16395*39*39*sizeof(double));

    /* Call f0_diag */
    xgc4py_f0_diag(f0_f_offset, f0_f_ndata, isp, 
                f0_f, f0_f_shap, f0_f_ndim,  /* f0_f data (input) */
                den, u_para, T_perp, T_para); /* (output) */

    /* Checking */
    double sum = 0.0;
    for (int i=0; i<16395*39*39; i++)
    {
        sum += den[i];
    }
    // Expected value: 493052035.5385171
    printf ("xgc4py_f0_diag: sum(den)= %f\n", sum);

    Py_Finalize();
    return 0;
}