#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int TestIGraph(void);

int main(void) {
  // Unit tests
  int resu=TestIGraph();
  if (resu) printf("IGraph test OK !\n");
  else printf("IGraph test failed !\n");
  return !resu;
} 

// some unit tests of the macromesh code
int TestIGraph(void)
{
  MacroMesh m;

  int deg[]={2, 2, 2};
  int raf[]={2, 2, 2};

  int test = true;
  
  // test gmsh file reading
  ReadMacroMesh(&m, "../test/testmacromesh.msh");
  //ReadMacroMesh(&m, "cubegros.msh");
  Detect2DMacroMesh(&m);
  BuildConnectivity(&m);
  CheckMacroMesh(&m, deg, raf);
  //PrintMacroMesh(&m);
  schnaps_real vit[3] = {1, 1, 0};
  BuildMacroMeshGraph(&m, vit, deg, raf);
  
  igraph_t* graph = &m.connect_graph;
  int nedges = igraph_ecount(graph);
  test = test && (nedges <= m.nbfaces);
  assert(test);

  int nverts = igraph_vcount(graph);
  test = test && (nverts == m.nbelems + 2);
  printf("nverts=%d\n",nverts);
  assert(nverts == m.nbelems + 2);
  igraph_bool_t is_dag;
  igraph_is_dag(graph, &is_dag);
  test = test && is_dag; 
  assert(is_dag);
  
  for(int eid = 0; eid < nedges; eid++){
    int ifa =  m.edge2face[eid];
    int ieL = m.face2elem[4 * ifa + 0];
    int ieR = m.face2elem[4 * ifa + 2];
    int iefrom, ieto;
    int dir = m.edge_dir[eid];
    igraph_edge(graph, eid, &iefrom, &ieto);
    printf("Edge %d from %d to %d. dir=%d ieL=%d ieR=%d\n",
	   eid, iefrom, ieto, dir, ieL,ieR);

    if (iefrom >= m.nbelems) iefrom = -1; 
    if (ieto >= m.nbelems) ieto = -1; 
    
    if (dir == 0) {
      test = test && (iefrom == ieL && ieto == ieR);
      assert(test);
    }
    if (dir == 1) {
      test = test && (iefrom == ieR && ieto == ieL);
      printf("pb iefrom=%d ieR=%d ieto=%d ieL=%d\n",iefrom,ieR,ieto,ieL);
      assert(iefrom == ieR && ieto == ieL);
    }
  }

  // print upwind and downwind neighbors
  igraph_vector_t found_vert;
  igraph_vector_init(&found_vert, 6);
  for(int nid = 0; nid < nverts; nid++){
    printf("nid=%d\n",nid);
    igraph_neighbors(graph, &found_vert, nid, IGRAPH_IN);
    int nup = igraph_vector_size(&found_vert);
    printf("Vert. %d has %d Upwind Vert.:",nid,nup);
    for(int ii=0; ii < nup; ii++) {
      int upvert = (int)VECTOR(found_vert)[ii];
      printf(" %d ", upvert);
      //assert(m.topo_order[nid] > m.topo_order[upvert]);
    }
    igraph_neighbors(graph, &found_vert, nid, IGRAPH_OUT);
    int ndown = igraph_vector_size(&found_vert);
    printf("\nVert. %d has %d Downwind Vert.:",nid,ndown);
    for(int ii=0; ii < ndown; ii++) {
      int downvert = (int)VECTOR(found_vert)[ii];
      printf(" %d ",downvert);
      //assert(m.topo_order[nid] < m.topo_order[downvert]);
    }
    printf("\n");
  }
  // drawing
  FILE *f = fopen("macromesh.dot","w");;
  igraph_write_graph_dot(graph, f);
  fclose(f);

  printf("For drawing the graph:\n");
  printf("dot -Tpdf macromesh.dot -o graph.pdf\n");


  return test;
}
