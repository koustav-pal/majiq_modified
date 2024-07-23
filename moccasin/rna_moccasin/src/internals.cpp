/**
 * internals.cpp
 */

#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/npy_3kcompat.h>
#include <numpy/ufuncobject.h>

#include <moccasin/InferParams.hpp>

// for functions with no extra data being passed in
static void *data[] = {nullptr};

static PyMethodDef ModuleMethods[] = {{nullptr, nullptr, 0, nullptr}};
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  // namespace for module
  "internals",
  // could maybe change first nullptr to create docstring. Apparently should be
  // made with PyDoc_STRVAR?
  nullptr, -1, ModuleMethods,
  nullptr, nullptr, nullptr, nullptr};

PyMODINIT_FUNC PyInit_internals(void) {
  PyObject *m = PyModule_Create(&moduledef);
  if (!m) { return nullptr; }
  // get instance dictionary for module, which we will attach gufuncs to
  PyObject *d = PyModule_GetDict(m);

  import_array();
  import_umath();

  using moccasin::InferParams;
  PyObject *infer_params = PyUFunc_FromFuncAndDataAndSignature(
      InferParams::funcs,
      data,
      InferParams::types,
      InferParams::ntypes,
      InferParams::nin,
      InferParams::nout,
      PyUFunc_None,
      InferParams::name,
      InferParams::doc,
      0,
      InferParams::signature);
  PyDict_SetItemString(d, InferParams::name, infer_params);
  Py_DECREF(infer_params);

  return m;
}
