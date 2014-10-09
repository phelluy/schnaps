#include "macromesh.h"
#include <stdio.h>
#include <assert.h>
#include "geometry.h"
#include "global.h"
#include "test.h"

int main(void) {
  
  // unit tests
  
  assert(TestMacroMesh());
  
  assert(TestGeometry());
	 
  assert(TestInterpolation());

  assert(TestModel());

   assert(TestField()); 

  assert(TestFieldDG());

   //assert(TestFieldRK2());

  printf("All tests OK !\n");

  return 0;
} 
