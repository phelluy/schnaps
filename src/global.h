// \mainpage
//\f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f$
// global variables for the SCHNAPS project
/*! \mainpage SCHNAPS: Solveur Conservatif Hyperbolique Non-linéaire Appliqué aux PlasmaS
 *
 * \section intro_sec Introduction
 *
Mode d'emploi

téléchargement (nécessite git)

git clone https://github.com/phelluy/SCHNAPS.git

se placer dans le dossier SCHNAPS

puis:

cmake .

make

ctest

une série de tests démarre. Il est conseillé de faire repasser ces
tests après toute modification du code.

Pour lancer schnaps:

générer le maillage:

gmsh disque.geo -3

./schnaps

Cet exemple consiste à résoudre l'équation de transport à
l'intérieur d'un cylindre. La visualisation  des résultats
utilise gmsh

gmsh dgvisu.msh

(puis tools -> options -> view[0] -> adapt visualization grid pour
afficher une image plus jolie...)

Adapter le fichier source "schnaps.c"
pour traiter des cas avec d'autres maillages. Des exemples de
maillages se trouvent dans le dossier geo.

Génération de la documentation (avec doxygen):

cd doc/
doxygen doxyschnaps
 *
 *
 */
#ifndef _GLOBAL_H
#define _GLOBAL_H

//! brief global variables and defs
#ifndef _OPENMP
// activate pthread if openmp is not here
//#define _WITH_PTHREAD
#endif

#define _WITH_OPENCL

extern int nplatform_cl;
extern int ndevice_cl;

char numflux_cl_name[1024];
char cl_buildoptions[1024];

#define __constant

#endif
