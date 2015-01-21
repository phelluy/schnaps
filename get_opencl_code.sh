# -*- mode: sh; -*-
#!/bin/sh
echo "#pragma start_opencl" > schnaps_temp.cl
echo "\n/*-----------------------------src from .h files-----------------------------*/\n" >> schnaps_temp.cl
echo "#pragma OPENCL EXTENSION cl_khr_fp64: enable" >> schnaps_temp.cl
echo "#pragma end_opencl" >> schnaps_temp.cl

cat src/*.h >> schnaps_temp.cl

echo "#pragma start_opencl" >> schnaps_temp.cl
echo "\n/*-----------------------------src from .c files-----------------------------*/\n" >> schnaps_temp.cl
echo "#pragma end_opencl" >> schnaps_temp.cl

cat src/*.c >> schnaps_temp.cl

echo "// -*- mode: c; -*-\n" > schnaps.cl

sed -n '/#pragma start_opencl/,/#pragma end_opencl/p' schnaps_temp.cl >> schnaps.cl
sed -i '/start_opencl/d' schnaps.cl
sed -i '/end_opencl/d' schnaps.cl

echo "\n/*-----------------------------src from .cl files----------------------------*/\n" >> schnaps.cl

cat src/field.cl >> schnaps.cl
#rm schnaps_temp.cl
return 0
