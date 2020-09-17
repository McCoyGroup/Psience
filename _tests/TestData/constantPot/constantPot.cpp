#include "Python.h"
#include <vector>
#include <string>

// All the needs to be changed is this

double constantPot(std::vector<std::vector<double> >, std::vector<std::string>) {
    return 52.0;
}

/******************** FIND AND REPLACE ON constantPot FROM HERE ON OUT ****************************/
// Everything else can be handled through find-replace
PyObject *_WrapFunction( void *ptr) {
    PyObject *link_cap = PyCapsule_New(ptr, "potential", NULL);
    if (link_cap == NULL) {
        PyErr_SetString(PyExc_TypeError, "couldn't create capsule object");
        return NULL;
    } else if (!PyCapsule_IsValid(link_cap, "potential")) {
        PyErr_SetString(PyExc_ValueError, "couldn't add pointer to invalid capsule object");
        Py_XDECREF(link_cap);
        return NULL;
    }
    return link_cap;
}

static PyMethodDef constantPotMethods[] = {
    {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION > 2

const char constantPot_doc[] = "exposes constantPot as a potential";
static struct PyModuleDef constantPotModule = {
    PyModuleDef_HEAD_INIT,
    "constantPot",   /* name of module */
    constantPot_doc, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
   constantPotMethods
};

PyMODINIT_FUNC PyInit_constantPot(void)
{
    PyObject * m = PyModule_Create(&constantPotModule);
    if ( m == NULL ) return NULL;
    PyModule_AddObject(m, "potential", _WrapFunction( (void *) *constantPot ));

    return m;
}
#else

PyMODINIT_FUNC initconstantPot(void)
{
    m = Py_InitModule("constantPot", constantPotMethods);
    if ( m == NULL ) return NULL;
    PyModule_AddObject(m, "potential", _WrapFunction( (void *) *constantPot ));
}

#endif
