#ifdef _WITH_PYTHONLIBS
#include "python_call.h"

int main(int argc, char *argv[])
{   
    char *ocl_buffer = NULL;
    Py_ssize_t ocl_buffer_size = 1;

    if (call_py_opencl_src_extract(&ocl_buffer,&ocl_buffer_size)) {
        return -1;
    } else {
        printf("%s\n",ocl_buffer);
        printf("Size = %d\n",ocl_buffer_size);
        assert(ocl_buffer[ocl_buffer_size]=='\0');
        return 0;
    }
}
#else
int main(int argc, char *argv[])
{
    return 0;
} 
#endif
