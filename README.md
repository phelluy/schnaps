# SCHNAPS

Solveur pour les lois de Conservation Hyperboliques Non-linéaires
Appliqué aux PlasmaS

## Mode d'emploi

### Téléchargement (nécessite git)

Accès développeur:
`git@gitlab.math.unistra.fr:tonus/schnaps.git`

Accès lecture seule:
`git clone https://gforge.inria.fr/git/schnaps/schnaps.git`

### Compilation

se placer dans le dossier schnaps

créer un dossier

	mkdir build

se placer dans ce dossier

	cd build

If StarPU isn't in the standard (system-wide) location, one can
specify the location with export STARPU_DIR=<...>

puis:

	cmake ..
	make
	ctest

Une série de tests démarre. Il est conseillé de faire repasser ces
tests après toute modification du code.

Si certains tests StarPU ne passent pas, désactiver les codelets
OpenCL et CUDA avant de lancer ctest:

	export STARPU_NOPENCL=0
	export STARPU_NCUDA=0

### Pour lancer schnaps:

#### générer le maillage:

	gmsh disque.geo -3

	./schnaps

Cet exemple consiste à résoudre l'équation de transport à
l'intérieur d'un cylindre.

(Edit: en ce moment cet exemple est désactivé) 

La visualisation  des résultats utilise gmsh

	gmsh dgvisu.msh

(puis tools -> options -> view[0] -> adapt visualization grid pour
afficher une image plus jolie...)

Adapter le fichier source `schnaps.c` pour traiter des cas avec
d'autres maillages. Des exemples de maillages se trouvent dans le
dossier geo.


#### SCHNAPS est sous licence CeCILL:

http://www.cecill.info/licences/Licence_CeCILL_V1.1-US.html

#### Génération de la documentation (avec doxygen):

	cd doc/
	doxygen doxyschnaps
	 *
	 *
	 */

#### Installation de StarPU (préalable à l'installation de schnaps)

Facultatif (permet de visualiser les traces starpu):
Télécharger FxT

	mkdir /usr/local/fxtdir
	cd FxT
	./bootstrap
	export FXTDIR=/usr/local/fxtdir
	./configure --prefix=$FXTDIR
	make -j4
	make install

Indispensable:
Télécharger StarPU a partir de la base svn (voir site StarPU)
brew 

	cd StarPu/
	mkdir build
	./autogen.sh
(installer les paquets brew manquants)
Si Fxt a été installé:

	../configure --with-fxt=$FXTDIR

sinon:

	../configure
	make -j4
	make check (facultatif)
	make install

Remarque: il y a une limitation de taille des kernels opencl dans
StarPU qui est corrigée avec ce patch:

	Index: `src/drivers/opencl/driver_opencl_utils.c

	--- src/drivers/opencl/driver_opencl_utils.c    (révision 16654)
	+++ src/drivers/opencl/driver_opencl_utils.c    (copie de travail)
	@@ -364,7 +364,7 @@
	     char located_file_name[1024];
	     char located_dir_name[1024];
	     char new_build_options[1024];
	-    char opencl_program_source[16384];
	+    char opencl_program_source[131072];

     // Do not try to load and compile the file if there is no devices
     nb_devices = starpu_opencl_worker_get_count();


Désactiver les codelets opencl et cuda (éventuellement) :

	export STARPU_NOPENCL=0 
	export STARPU_NCUDA=0 

Pour visualiser les traces, télécharger et compiler VITE a partir de
la base svn (voir site http://vite.gforge.inria.fr/) Installer
graphviz pour visualiser les graphes de tâche.


