#include "python_call.h"

/* Example for casting a int C array to PyObject */
PyObject *make_integer_list(long int *array, unsigned long size) {
    PyObject *py_list = PyList_New(size);
    for (unsigned int i = 0; i != size; ++i) {
        PyList_SET_ITEM(py_list, i, PyLong_FromLong(array[i]));
    }
    return py_list;
}


/* Example for casting a double C array to PyObject */
PyObject *make_double_list(double *array, unsigned long size) {
    PyObject *py_list = PyList_New(size);
    for (unsigned int i = 0; i != size; ++i) {
        PyList_SET_ITEM(py_list, i, PyLong_FromDouble(array[i]));
    }
    return py_list;
}

/* Call function extract_ocl from ocl_extract module*/

int call_py_opencl_src_extract(
        char **ocl_buffer,
        Py_ssize_t *ocl_buffer_size)
{
    PyObject *pName = NULL;
    PyObject *pModule = NULL;
    PyObject *pFunc = NULL;
    PyObject *pValue = NULL;

    char *python_file;
    char *python_function;
    // const char *ocl_extracted_buffer = NULL;
    
    python_file = "ocl_extract";
    python_function = "extract_ocl";

    Py_Initialize();
    /* Hack to include another python path..*/
    PyRun_SimpleString("import sys\nsys.path.append(\"../python\")\n");
#ifdef IS_PY3
    pName = PyUnicode_DecodeFSDefault(python_file);
#else 
    pName = PyUnicode_DecodeUTF8(python_file, strlen(python_file) , NULL);
#endif
    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        
        pFunc = PyObject_GetAttrString(pModule, python_function);

        if (pFunc && PyCallable_Check(pFunc)) {
            /* Python function call */
            pValue = PyObject_CallObject(pFunc, NULL);
            
            if (pValue != NULL) {
                /* Decoding and construct C char OpenCL buffer */
#ifdef IS_PY3
                const char *res = PyUnicode_AsUTF8AndSize(pValue, ocl_buffer_size);
                *ocl_buffer = (char*)calloc(*ocl_buffer_size+1, sizeof(char)); 
                strncpy(*ocl_buffer ,res, *ocl_buffer_size+1);
#else
                // *ocl_buffer = PyString_AsString(pValue);
                int ret;
                char *res = NULL;
                ret = PyString_AsStringAndSize( pValue, 
                                                &res, 
                                                ocl_buffer_size);
                                                
                *ocl_buffer = (char*)calloc(*ocl_buffer_size, sizeof(char)); 
                strncpy(*ocl_buffer, res, *ocl_buffer_size);
                                        
#endif
                if (ocl_buffer == NULL) {
                    if (PyErr_Occurred()) PyErr_Print();
                    fprintf(stderr, "[Error] failed to decode\n");  
                    return 1;
                }
                Py_DECREF(pValue);
            } else {
                if (PyErr_Occurred()) PyErr_Print();
                fprintf(stderr, "[Error] failed function call\n");
                return 1;
            }
            Py_DECREF(pFunc);
        } else {
            if (PyErr_Occurred()) PyErr_Print();
            fprintf(stderr, 
                    "[Error] failed to find function \"%s\"\n", 
                    python_function);
            return 1;
        }
        
        /* Free memory */
        Py_DECREF(pModule);
    } else {
        if (PyErr_Occurred()) PyErr_Print();
        fprintf(stderr, "[Error] failed to load module \"%s\"\n", python_file);
        return 1;
    }
    /* Kill interpreter*/
    Py_Finalize();
    
    return 0;
}