#include "macromesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>


void ReadMacroMesh(MacroMesh* m,char* filename){

  FILE* f=NULL;

  char* line=NULL;
  size_t linesize=0;
  size_t ret;

  f=fopen(filename,"r");
  assert(f != NULL);

  do {
    getline(&line,&linesize,f);
  }
  while(strcmp(line,"$Nodes\n") != 0);   

  // read the nodes data
  getline(&line,&linesize,f);
  m->nbnodes=atoi(line);
  printf("nbnodes=%d\n",m->nbnodes);

  m->node=malloc(3 * m->nbnodes * sizeof(double));
  assert(m->node);

  for(int i=0;i<m->nbnodes;i++){    
    getdelim(&line,&linesize,(int) ' ',f); // node number
    getdelim(&line,&linesize,(int) ' ',f); // x
    m->node[3*i+0]=atof(line);
    getdelim(&line,&linesize,(int) ' ',f); // y
    m->node[3*i+1]=atof(line);
    getline(&line,&linesize,f); // z (end of the line)
    m->node[3*i+2]=atof(line);
    /* printf("Node %d x=%f y=%f z=%f\n",i, */
    /* 	   m->node[3*i+0], */
    /* 	   m->node[3*i+1], */
    /* 	   m->node[3*i+2]); */
  }    
  getline(&line,&linesize,f);
  //printf("%s",line);
  // check that we have reached the end of nodes
  assert(strcmp(line,"$EndNodes\n") == 0);
  
  // Now read all the elements of the mesh
  do {
    getline(&line,&linesize,f);
  }
  while(strcmp(line,"$Elements\n") != 0);   

  // size of the gmsh elems list 
  getline(&line,&linesize,f);
  int nball=atoi(line);
  // allocate to a too big size
  m->elem2node=malloc(20 * sizeof(int) * nball);
  assert(m->elem2node);

  // now count only the H20 elems (code=17)
  m->nbelems=0;
  int countnode=0;
  for(int i=0;i<nball;i++){
    getdelim(&line,&linesize,(int) ' ',f); // elem number
    getdelim(&line,&linesize,(int) ' ',f); // elem type
    int elemtype=atoi(line);
    if (elemtype != 17) {
      getline(&line,&linesize,f);
    }
    else {
      m->nbelems++;
      getdelim(&line,&linesize,(int) ' ',f); //useless code
      getdelim(&line,&linesize,(int) ' ',f); //useless code
      getdelim(&line,&linesize,(int) ' ',f); //useless code
      for(int j=0;j<19;j++){
	getdelim(&line,&linesize,(int) ' ',f);
	//printf("%d ",atoi(line));
	m->elem2node[countnode]=atoi(line);
	countnode++;
      }
      getline(&line,&linesize,f);
      //printf("%d\n",atoi(line));
      m->elem2node[countnode]=atoi(line)-1;
      countnode++;
    }
  }	
  getline(&line,&linesize,f);
  //printf("%s",line);
  // check that we have reached the end of nodes
  assert(strcmp(line,"$EndElements\n") == 0);
  printf("nbelems=%d\n",m->nbelems);
  realloc(m->elem2node,20 * sizeof(int) * m->nbelems);
  assert(m->elem2node);

}

// display macromesh dat on standard output
void PrintMacroMesh(MacroMesh* m){
  printf("nbnodes=%d\n",m->nbnodes);
  for(int i=0;i<m->nbnodes;i++){
    printf("node %d x=%f y=%f z=%f\n",i,
	   m->node[3*i+0],
	   m->node[3*i+1],
	   m->node[3*i+2]);
  }
  printf("nbelems=%d\n",m->nbelems);
  for(int i=0;i<m->nbelems;i++){
    printf("elem %d -> ",i);
    for(int j=0;j<20;j++){
      printf("%d ",m->elem2node[20*i+j]);
    }
    printf("\n");
  }
}
