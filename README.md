SCHNAPS
=======

Solver for Conservative Hypebolic Non-linear systems Applied to PlasmaS

Solveur pour les lois de Conservation Hyperboliques Non-linéaires
Appliqué aux PlasmaS


Downloads:

Develloper access:
git clone git+ssh://<gforge_account_name>\@scm.gforge.inria.fr//gitroot/schnaps/schnaps.git

Read-only access:
git clone https://gforge.inria.fr/git/schnaps/schnaps.git


Compilation:

From the schnaps folder, first

mkdir build

cd build

cmake ..

make

gmsh ../geo/disque.geo -3

ctest


By default, the OpenCL code runs on platform 0, device 0.  This
can be changed in schnaps via command-line argument.  The unit tests
run on the default platform/device, which can be changed via

cmake -D_CL_PLATFORM=1 -D_CL_DEVICE=1  ..

By default, computations are done in double-precision.  This can be
changed via

cmake -DSINGLE_PRECISION:BOOL=ON


The main schnaps program is schnaps.  To generate .msh files from .geo
files, one run

gmsh <file>.geo -3

The resulting .msh file can be passed to schnaps via command-line
arguments (see ./schnaps -h for a list of options).

The default arguments for schnaps results in a simulation of 2D
transport in a disc.  The results can be viewed via

gmsh dgvisu.msh

The default gmsh visualizer is quite coarse, but the visualization can
be improved by the setting "Adapt visuzliation grid" in
tools -> options -> View [0].


There is some basic documentation with doxygen, which can be generated via

cd doc/
doxygen doxyschnaps



SCHNAPS is under the CeCILL license:
http://www.cecill.info/licences/Licence_CeCILL_V1.1-US.html
