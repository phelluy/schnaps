#!/bin/sh
cat src/*.h > schnaps_temp.cl
cat src/*.c >> schnaps_temp.cl
sed -n '/#pragma start_opencl/,/#pragma end_opencl/p' schnaps_temp.cl > schnaps.cl
sed -i '/start_opencl/d' schnaps.cl
sed -i '/end_opencl/d' schnaps.cl
cat src/field.cl >> schnaps.cl
return 0
