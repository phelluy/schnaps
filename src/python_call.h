/* Python 3 embedding routines for OpenCl code generation.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 and
 * only version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not see <http://www.gnu.org/licenses/>.
 */

#ifndef _PYTHON_CALL_
#define _PYTHON_CALL_

#include <Python.h>
#include <stdlib.h>

#if PY_MAJOR_VERSION >= 3
#define IS_PY3
#endif


PyObject *make_integer_list(
        long int *array, 
        unsigned long size);
        
PyObject *make_double_list(
        double *array,
        unsigned long size);

int call_py_opencl_src_extract(
        char **ocl_buffer, 
        Py_ssize_t *ocl_buffer_size);

#endif