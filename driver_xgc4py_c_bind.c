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

    double data[3] = {1.,2.,3.};
    xgc4py_test_print(data, 3);
    Py_Finalize();
    return 0;
}