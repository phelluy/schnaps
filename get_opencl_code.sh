#!/bin/sh

rm -f schnaps_temp.cl
touch schnaps_temp.cl

echo "#pragma start_opencl" >> schnaps_temp.cl
echo "#ifdef cl_khr_fp64" >> schnaps_temp.cl
echo "#pragma OPENCL EXTENSION cl_khr_fp64: enable" >> schnaps_temp.cl
echo "#else" >> schnaps_temp.cl
echo "#error" >> schnaps_temp.cl
echo "#endif" >> schnaps_temp.cl
echo "#pragma end_opencl" >> schnaps_temp.cl

# NB: echo "\n/* does not work on the mesocentre.

cat src/*.h >> schnaps_temp.cl
cat src/*.c >> schnaps_temp.cl

# Copy everything between "#pragma start_opencl" and "#pragma end_opencl"
sed -n '/#pragma start_opencl/,/#pragma end_opencl/p' schnaps_temp.cl > schnaps.cl

## Remove undefined pragmas
cat schnaps.cl | sed 's/\#pragma\ end_opencl//' > schnaps0.cl
mv schnaps0.cl schnaps.cl

cat schnaps.cl | sed 's/\#pragma\ start_opencl//' > schnaps0.cl
mv schnaps0.cl schnaps.cl

cat src/field.cl >> schnaps.cl

rm -f schnaps_temp.cl

exit 0
