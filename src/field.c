#include "field.h"


// param[0]=

// param[0] = M
// param[1] = deg x
// param[2] = deg y
// param[3] = deg z
// param[4] = raf x
// param[5] = raf y
// param[6] = raf z
int GenericVarindex(int* param, int elem, int ipg, int iv){

  int npg= (param[1]+1)*(param[2]+1)*(param[3]+1);

  return iv + param[0] * ( ipg + npg * elem);

}
