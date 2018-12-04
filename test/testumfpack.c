#include "solverumfpack.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "test.h"

int TestUmfPack(void);

int main(void) {
  // unit tests
  int resu=TestUmfPack();
  if (resu) printf("UmfPack test OK !\n");
  else printf("UmfPack test failed !\n");
  return !resu;
} 

int TestUmfPack(void){
  int test = true;
  // reference element
  smalltestumfpack();
  return test;
}
