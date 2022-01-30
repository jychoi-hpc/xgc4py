#include <Python.h>
#include "xgc4py_c_bind.h"

int
main() 
{
    PyImport_AppendInittab("xgc4py_c_bind", PyInit_xgc4py_c_bind);
    Py_Initialize();
    PyImport_ImportModule("xgc4py_c_bind");
    call_hello(77);
    Py_Finalize();
    return 0;
}