#!/bin/sh
#exit
rm -f schnaps.cl
touch schnaps.cl

echo "#ifdef cl_khr_fp64" >> schnaps.cl
echo "#pragma OPENCL EXTENSION cl_khr_fp64: enable" >> schnaps.cl
echo "#else" >> schnaps.cl
echo "//#error" >> schnaps.cl
echo "#endif" >> schnaps.cl
echo "#define M_PI 3.1415926535897932385" >> schnaps.cl


echo >> schnaps.cl

echo "#define NULL 0" >> schnaps.cl

# NB: echo "\n/* does not work on the mesocentre.

# Copy everything between "#pragma start_opencl" and "#pragma end_opencl"
sed -n '/#pragma start_opencl/,/#pragma end_opencl/p' ../src/*.h  >> schnaps.cl
sed -n '/#pragma start_opencl/,/#pragma end_opencl/p' ../src/*.c  >> schnaps.cl


# Remove undefined pragmas
cat schnaps.cl | sed 's/\#pragma\ end_opencl//' > schnaps0.cl
mv schnaps0.cl schnaps.cl
cat schnaps.cl | sed 's/\#pragma\ start_opencl//' > schnaps0.cl
mv schnaps0.cl schnaps.cl

cat ../src/field.cl >> schnaps.cl

rm -f schnaps0.cl

exit 0
	 #  3 - test/testfieldrk2_cl (OTHER_FAULT)
	 #  8 - test/testcharge_ocl (SEGFAULT)
	 # 24 - test/testfielddg_spu (OTHER_FAULT)
	 # 29 - test/testmeq2_cl (OTHER_FAULT)
	 # 44 - test/testfieldrk4_cl (OTHER_FAULT)
	 # 51 - test/testlinearsolver (Failed)
