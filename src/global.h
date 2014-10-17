// \mainpage 
//\f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f$
// global variables for the SCHNAPS project
/*! \mainpage SCHNAPS: Solveur Conservatif Hyperbolique Non-linéaire Appliqué aux PlasmaS
 *
 * \section intro_sec Introduction
 *
 * SCHNAPS est un solveur générique Galerkin Discontinu pour les lois de conservation
 *
 * \section install_sec Installation
 *
 * git clone https://github.com/phelluy/SCHNAPS.git

se placer dans le dossier SCHNAPS

puis:

cmake .

make

ctest

une série de tests démarre. Il est conseillé de faire repasser ces
tests après toute modification du code.

Pour lancer schnaps:

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
Pour générer un maillage à partir du fichier disque.geo par exemple,
faire:

gsmh disque.geo -3

Cette commande génère un fichier disque.msh lisible par SCHNAPS.

Génération de la documentation (avec doxygen):

cd doc/

doxygen

doxyschnaps

 *  
 * 
 */
#ifndef _GLOBAL_H
#define _GLOBAL_H

// precision: replace by float for single precision
//#define double double
//#define double float

// suffix for floating point constant
// nothing for double precision
// f for single precision
#define _FP
//#define _FP f

//#define NULL ((void*) 0)

#define _DEGX 3
#define _DEGY 3
#define _DEGZ 3
#define _RAFX 1
#define _RAFY 1
#define _RAFZ 1



#endif
