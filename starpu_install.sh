#!/bin/sh

# Locally install FXT and StarPU on irma-atlas
# Works for me (Michel Massaro)

VERT="\\033[1;32m"
NORMAL="\\033[0;39m"
ROUGE="\\033[1;31m"

INSTALL_FXT_DIR=$HOME
INSTALL_SPU_DIR=$HOME


echo "$ROUGE" "Install FXT in '$INSTALL_FXT_DIR/fxt-0.2.11'" "$NORMAL"

cd $INSTALL_FXT_DIR

wget http://download.savannah.gnu.org/releases/fkt/fxt-0.2.11.tar.gz
tar -xvzf fxt-0.2.11.tar.gz
rm fxt-0.2.11.tar.gz
cd fxt-0.2.11/
./configure --prefix=$INSTALL_FXT_DIR/fxt-0.2.11/
make -j 8
make install

echo "$ROUGE" "Install StarPU in '$INSTALL_SPU_DIR/StarPU'" "$NORMAL"

cd $INSTALL_SPU_DIR

svn checkout svn://scm.gforge.inria.fr/svn/starpu/trunk StarPU
cd StarPU
mkdir build



# Now the fix is implemented

## Fix bug before the installation
#echo "$VERT" "Apply the following modification ..." "$NORMAL"
#
#sed -i -e '/_opencl_compile_or_load_opencl_from_string/,$ {s/static char buffer/\/\/static char buffer/; ta}; b; :a; N; P; s/^[^\n]*\n//; ba' src/drivers/opencl/driver_opencl_utils.c
#sed -i -e '/_opencl_compile_or_load_opencl_from_string/,$ {s/clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, \&len);/clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, \&len);\n                        char *buf2 = malloc(2 * (len + 1));\n                        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, len, buf2, NULL);/; ta}; b; :a; N; P; s/^[^\n]*\n//; ba' src/drivers/opencl/driver_opencl_utils.c
#sed -i-e '/_opencl_compile_or_load_opencl_from_string/,$ {s/_STARPU_DISP("Compilation output\\n%s\\n", buffer);/_STARPU_DISP("Compilation output\\n%s\\n", buf2);\n                        free(buf2);/; ta}; b; :a; N; P; s/^[^\n]*\n//; ba' src/drivers/opencl/driver_opencl_utils.c
#
#sed -i -e '/_opencl_compile_or_load_opencl_from_file/,$ {s/char opencl_program_source\[16384\]/char opencl_program_source\[200000\]/; ta}; b; :a; N; P; s/^[^\n]*\n//; ba' src/drivers/opencl/driver_opencl_utils.c
#
#echo "$VERT" "Done" "$NORMAL"




# Enable subdir_objects for automake 1.14.*
sed -i 's/foreign/foreign subdir-objects/g' configure.ac

./autogen.sh
cd build
../configure --with-fxt=$INSTALL_FXT_DIR/fxt-0.2.11/ --prefix=$INSTALL_SPU_DIR/StarPU/ --disable-build-doc

# Disable example compilation because there is a bug
sed -i 's/\$(am__append_3) //g' Makefile
make -j8
#make check
make install

cd $HOME

echo ""
echo "$VERT" "To finish, you need to add some linking" "$NORMAL"
echo "$VERT" "In my .bashrc file, I have added these line" "$NORMAL"

echo "     export PKG_CONFIG_PATH=$INSTALL_SPU_DIR/StarPU/lib/pkgconfig"
echo "     PKG_CONFIG_PATH=\$PKG_CONFIG_PATH:$INSTALL_SPU_DIR/StarPU"
echo "     LD_LIBRARY_PATH=$INSTALL_SPU_DIR/StarPU/lib:$LD_LIBRARY_PATH"
echo "     export STARPU_INCLUDE_DIRS=$INSTALL_SPU_DIR/StarPU/include"
echo "     export STARPU_LIBRARY_DIRS=$INSTALL_SPU_DIR/StarPU/lib"

echo ""
echo "$VERT" "It is also necessary to specify the directory for compiling Schnaps" "$NORMAL"
echo "$VERT" "I use the shortcut" "$NORMAL"
echo "     alias cmakespu='cmake .. -DSTARPU=$INSTALL_SPU_DIR//StarPU'"



echo ""
echo ""
echo "$VERT" "I can do it for you [y/n]" "$NORMAL"
read resp
if [ $resp = "y" ]
then
    echo "" >> ~/.bashrc
    echo "" >> ~/.bashrc
    echo "export PKG_CONFIG_PATH=$INSTALL_SPU_DIR/StarPU/lib/pkgconfig" >> ~/.bashrc
    echo "PKG_CONFIG_PATH=\$PKG_CONFIG_PATH:$INSTALL_SPU_DIR/StarPU" >> ~/.bashrc
    echo "LD_LIBRARY_PATH=$INSTALL_SPU_DIR/StarPU/lib:$LD_LIBRARY_PATH" >> ~/.bashrc
    echo "export STARPU_INCLUDE_DIRS=$INSTALL_SPU_DIR/StarPU/include" >> ~/.bashrc
    echo "export STARPU_LIBRARY_DIRS=$INSTALL_SPU_DIR/StarPU/lib" >> ~/.bashrc
    echo "alias cmakespu='cmake .. -DSTARPU=$INSTALL_SPU_DIR//StarPU'" >> ~/.bashrc

    source ~/.bashrc
fi
