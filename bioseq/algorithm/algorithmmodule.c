#include "algorithm.h"
#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject *
algorithm_NeedlemanWunsch(PyObject *self, PyObject *args)
{
    char *query;
    char *subject;
    float match;
    float mismatch;
    float gap_open;
    float gap_extend;

    if (!PyArg_ParseTuple(args, "ssffff", &query, &subject, &match, &mismatch, &gap_open, &gap_extend))
        Py_RETURN_NONE;

    int size = strlen(query) + strlen(subject);
    char *align_query = malloc(size);
    char *align_subject = malloc(size);
    float score;
    NeedlemanWunsch(query, subject, align_query, align_subject, &score, match, mismatch, gap_open, gap_extend);

    PyObject *result = PyTuple_Pack(3, PyUnicode_FromString(align_query), PyUnicode_FromString(align_subject), PyFloat_FromDouble(score));

    free(align_query);
    free(align_subject);
    return result;
}

static PyObject *
algorithm_SmithWaterman(PyObject *self, PyObject *args)
{
    char *query;
    char *subject;
    float match;
    float mismatch;
    float gap_open;
    float gap_extend;

    if (!PyArg_ParseTuple(args, "ssffff", &query, &subject, &match, &mismatch, &gap_open, &gap_extend))
        Py_RETURN_NONE;

    int size = strlen(query) + strlen(subject);
    char *align_query = malloc(size);
    char *align_subject = malloc(size);
    float score;
    SmithWaterman(query, subject, align_query, align_subject, &score, match, mismatch, gap_open, gap_extend);

    PyObject *result = PyTuple_Pack(3, PyUnicode_FromString(align_query), PyUnicode_FromString(align_subject), PyFloat_FromDouble(score));

    free(align_query);
    free(align_subject);
    return result;
}

static PyMethodDef AlgorithmMethods[] = {
    {"NeedlemanWunsch", algorithm_NeedlemanWunsch, METH_VARARGS, "algorithm NeedlemanWunsch."},
    {"SmithWaterman", algorithm_SmithWaterman, METH_VARARGS, "algorithm SmithWaterman."},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef algorithmmodule = {
    PyModuleDef_HEAD_INIT,
    "algorithm",
    "algorithm module writed by c program language.",
    -1,
    AlgorithmMethods,
};

PyMODINIT_FUNC
PyInit_algorithm(void)
{
    return PyModule_Create(&algorithmmodule);
}