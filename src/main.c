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

  return 0;
} 
