// \mainpage
//\f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f$
// global variables for the SCHNAPS project
/*! \mainpage SCHNAPS: Solveur Conservatif Hyperbolique Non-linéaire Appliqué aux PlasmaS
 *
 * \section intro_sec Introduction
 *

SCHNAPS est sous licence CeCILL:

http://www.cecill.info/licences/Licence_CeCILL_V1.1-US.html

Mode d'emploi

téléchargement (nécessite git)

Accès développeur:
git clone git+ssh://<gforge_account_name>\@scm.gforge.inria.fr//gitroot/schnaps/schnaps.git

Accès lecture seule:
git clone https://gforge.inria.fr/git/schnaps/schnaps.git

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

#ifdef _WITH_OPENCL

extern int nplatform_cl;
extern int ndevice_cl;

char numflux_cl_name[1024];
char cl_buildoptions[1024];


#endif //_WITH_OPENCL
#define __constant
#define __local
#endif

#define _DOUBLE_PRECISION
#ifdef _DOUBLE_PRECISION
#define real double
#else
#define real float
#endif

#define __constant
