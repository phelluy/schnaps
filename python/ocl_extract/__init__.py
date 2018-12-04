import logging
import os
# import time
import sys


log_file_name = "extract_ocl.log"
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

handler = logging.FileHandler(log_file_name, mode = 'w')
handler.setLevel(logging.INFO)
# format = logging.Formatter("%(asctime)s - %(funcName)s - %(levelname)s - "
                           # + "%(message)s")
format = logging.Formatter("%(funcName)s - %(levelname)s - "
                           + "%(message)s")
handler.setFormatter(format)
logger.addHandler(handler)


def insert_opencl_header():
    ocl_header =  "#ifdef cl_khr_fp64\n"
    ocl_header += "#pragma OPENCL EXTENSION cl_khr_fp64: enable\n"
    ocl_header += "#else\n"
    ocl_header += "//#error\n"
    ocl_header += "#endif\n"
    ocl_header += "#define M_PI 3.1415926535897932385\n"
    ocl_header += "\n"
    ocl_header += "#define NULL 0\n"
    ocl_header += "\n"
    return ocl_header




def extract_opencl_pragma_file(file_name):
    '''
    
    Extract source code between "start_flag" and "end_flag".
    Return a string containing the extracted source code
    
    '''
    start_flag = "#pragma start_opencl"
    end_flag = "#pragma end_opencl"
    extract_str = ""
    if (version_major > 2):
        f = open(file_name, 'r', encoding = 'utf8')
    else :
        f = open(file_name, 'r')
    start_idx = None
    end_idx = None
    
    try:
        #Convert file to python list
        src_list = f.read().splitlines()
    
        #To avoid right useless space. 
        # src_list = [line.rstrip() for line in src_list]
    
        #Get index of lines with ocl pragma
        start_idx = [k for k, p in enumerate(src_list) if p == start_flag]
        end_idx = [k for k, p in enumerate(src_list) if p == end_flag]
            
        if(start_idx != []):
            logger.info(file_name)
            logger.info("    #pragma start_opencl : [" 
                        + ",".join(str(i) for i in start_idx) 
                        + "] len = " + str(len(start_idx)))
                        
            logger.info("    #pragma end_opencl   : [" 
                        + ",".join(str(i) for i in end_idx) 
                        + "] len = " + str(len(end_idx)))
                        
            #To avoid open but not closed pragma
            #TODO assert end_idx[k] > start_idx[k] for all k
            len_start = len(start_idx)
            len_end = len(end_idx)
            len_min = min(len_start,len_end)
        
  
            extract_list = [src_list[start_idx[i]+1:end_idx[i]] 
                            for i in range(0,len_min)]
  
                            
            extract_str = "\n\n".join( "\n".join(j for j in i) 
                            for i in extract_list)
        
        
    finally:
        if f is not None:
            f.close()
        
    return extract_str
 

 
def extract_opencl_pragma_folder(folder_name):
    '''
    
    Extract OpenCl source code, from all .c and .h files in "folder_name".
    OpenCL code is sorted such as : first .h file content followed 
    by .c file content.
    Return a string containing the extracted source code
    
    '''
    h_file = []
    c_file = []
  
    extract_str = ""
    
    #Get all .h files and .c files list
    for subdir, dirs, files in os.walk(folder_name):
        for file in files :
            ext = os.path.splitext(file)[-1].lower()
            file_name = os.path.join(subdir, file)
            if ext == ".h":
                h_file.append(file_name)
            if ext == ".c":
                c_file.append(file_name)
 
                
    
    h_file.sort()
    c_file.sort()
    
    #Extraction of OpenCL code from .h files
    for file_name in h_file:
        extracted_file = extract_opencl_pragma_file(file_name)
        if extracted_file != "":
            extract_str += extracted_file + "\n\n"
    
    #Extraction of OpenCL code from .c files
    for file_name in c_file:
        extracted_file = extract_opencl_pragma_file(file_name)
        if extracted_file != "":
            extract_str += extracted_file + "\n\n"
    
    return extract_str



def extract_opencl_file(file_name):
    try:
        if (version_major > 2):
            f = open(file_name, 'r', encoding = 'utf8')
        else :
            f = open(file_name, 'r')
        src_list = f.read()
    finally:
        f.close()
        return src_list
    

    
    
def extract_ocl():

    global version_major
    version_major = sys.version_info.major
    
    #Call the parser
    folder_name = "../src"
    field_file_name = "../src/field.cl"

    # start = time.time() 
    
    open_cl_file = insert_opencl_header()
    open_cl_file += extract_opencl_pragma_folder(folder_name)
    open_cl_file += extract_opencl_file(field_file_name)
    
    
    if (version_major > 2):
        ocl_file_ptr = open("schnaps.cl", "w", encoding = 'utf8')
    else:
        ocl_file_ptr = open("schnaps.cl", "w")
        
    ocl_file_ptr.write(open_cl_file)
    ocl_file_ptr.close()
    return open_cl_file
    # print("Execution time : ", time.time() - start)
