#include "implicit.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "waterwave2d.h"

schnaps_real vit[3] = {-1.4, -0.7, 0};

void KLU_Upwind_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux);
  void TestSteady_KLU_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w);
void TestSteady_KLU_InitData(schnaps_real *x, schnaps_real *w);
void TestSteady_KLU_Source(const schnaps_real *xy, const schnaps_real t, const schnaps_real *w, schnaps_real *S);
void KLU_Steady_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
			      schnaps_real *flux);

int test_klu_basic_orig(void);
int test_klu_basic(void);
int test_KLU_block_ut(int n, int block_size);
void test_KLU_linsolv(int n,int block_size);
int Test_KLU_Steady(void);

int main(void) {
/*  test_KLU_linsolv(8,2);*/
/*  assert(1==2);*/
  //test_klu_basic_orig();
  //
  //int t0= test_klu_basic();
  //
  //assert(false);
/*  test_KLU_block_ut(8,3);*/
/*  assert(2==3);*/
  //
/*  int n=5;*/
/*  KLU klumat;*/
/*  //*/
/*  InitKLU(&klumat, n);*/
/*  for(int i = 0; i < n; i++){*/
/*    int j = i + 1;*/
/*    if (j < n) SwitchOnKLU(&klumat, i, j);*/
/*    SwitchOnKLU(&klumat, i, i);*/
/*  }*/
/*  SwitchOnKLU(&klumat, 3, 4);*/

/*  printf("construct CSR struct\n");*/
/*  AllocateKLU(&klumat);*/
/*  printf(" KLU found %i blocks , %i offdiagblocks , largest block size %i\n",*/
/*    klumat.symbolic->nblocks,klumat.symbolic->nzoff,klumat.symbolic->maxblock);*/
/*  //*/
/*  printf("copy CSR struct\n");*/
/*  AllocateCopyKLU(&klumat);*/
/*  //*/
/*  assert(1==2);*/
  // unit test
  int resu1 = 0;
  resu1=Test_KLU_Steady();
  if (resu1 >=1) printf("KLU  test OK !\n");
  else printf("KLU test failed !\n");

  return !resu1;
} 
// sample code from klu documentation
int test_klu_basic_orig(void){
  int n = 5 ;
  int Ap [ ] = {0, 2, 5, 9, 10, 12} ;
  int Ai [ ] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
  double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
  double b [ ] = {8., 45., -3., 3., 19.} ;
  klu_symbolic *Symbolic ;
  klu_numeric *Numeric ;
  klu_common Common ;
  int i ;
  klu_defaults (&Common) ;
  Symbolic = klu_analyze (n, Ap, Ai, &Common) ;
  Numeric = klu_factor (Ap, Ai, Ax, Symbolic, &Common) ;
  klu_solve (Symbolic, Numeric, 5, 1, b, &Common) ;
  klu_free_symbolic (&Symbolic, &Common) ;
  klu_free_numeric (&Numeric, &Common) ;
  for (i = 0 ; i < n ; i++) printf ("x [%d] = %g\n", i, b [i]) ;
  
  return 1;
}
int test_klu_basic(void){
  int n = 5 ;
  int Ap [ ] = {0, 2, 5, 9, 10, 12} ;
  int Ai [ ] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
  schnaps_real Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
  schnaps_real xref[ ] = {1.,2.,3.,4.,5.};
  schnaps_real y [ ] = {8., 45., -3., 3., 19.} ;
  schnaps_real b [ ] = {8., 45., -3., 3., 19.} ;
  //
  KLU klumat;
  InitKLU(&klumat, n);
  //
  for (int j=0;j < n;j++){
    for (int iloc=Ap[j];iloc<Ap[j+1];iloc++){
      SwitchOnKLU(&klumat,Ai[iloc],j);
    }
  }
  //
  AllocateKLU(&klumat);
  DisplayKLU(&klumat);
  //
  for (int j=0;j < n;j++){
    for (int iloc=Ap[j];iloc<Ap[j+1];iloc++){
      SetKLU(&klumat,Ai[iloc],j,Ax[iloc]);
    }
  }
  //
  printf("Testing matrix/vector prod \n");
  MatVectKLU(&klumat,xref,y);
  for (int i=0;i < n;i++){
    printf(" y[%i] : %f  \t b[%i] : %f\n",i,y[i],i,b[i]);
    assert(y[i]==b[i]);
  }
  //
  printf("Numerical factorization\n");
  FactoKLU(&klumat);
  //
  SolveKLU(&klumat,b,y);
  //
  for (int i=0;i< n;i++){
    printf("xref[%i] : %f \t sol[%i] : %f \n",i,xref[i],i,y[i]);
  }
  //
  FreeKLU(&klumat);
  return 1;
}
//
void test_KLU_linsolv(int n,int block_size){
  //
  int seed=4043498;
  srand(seed);
  //
  int *PL = calloc(n,sizeof(int)); // line permutation vector
  int *PC = calloc(n,sizeof(int)); // column permutation vector
  int *PLinv = calloc(n,sizeof(int)); // line permutation vector (inverse)
  int *PCinv = calloc(n,sizeof(int)); // column permutation vector (inverse)

  // set permutation to identity
  for (int i=0;i< n;i++){
    PL[i]=i;
    PC[i]=i;
  }
  // generate random perturbation
  for (int i=0;i<n;i++){
    int r=rand()%n;
    int tmp=PL[r];
    PL[r]=PL[i];
    PL[i]=tmp;
  }  
  for (int i=0;i<n;i++){
    int r=rand()%n;
    int tmp=PC[r];
    PC[r]=PC[i];
    PC[i]=tmp;
  }
  for (int i=0;i< n;i++){
    PLinv[PL[i]]=i;
    PCinv[PC[i]]=i;
  }
  //
  LinearSolver lsol;
  MatrixStorage ms = KLU_CSR;
  Solver sv = LU;
  printf(" LinearSolver/KLU test with block triangular matrix\n");
  InitLinearSolver(&lsol,n,&ms,&sv);
  //
  for (int i=0;i < n; i++){
    for (int j=0; j< n; j++){
      if (i/block_size <= j/block_size){
        IsNonZero(&lsol,PL[i],PC[j]);
      }
    }
  }
  //
  AllocateLinearSolver(&lsol);
  //
  for (int i=0;i < n; i++){
    for (int j=0; j< n; j++){
      if (i/block_size <= j/block_size){
        schnaps_real tmp = ((schnaps_real) rand())/((schnaps_real) RAND_MAX); 
        SetLinearSolver(&lsol,PL[i],PC[j],tmp);
      }
    }
  }
  // 
  // generate solution 
  schnaps_real *xref = (schnaps_real*) calloc(n,sizeof(schnaps_real));
  schnaps_real *x = (schnaps_real*) calloc(n,sizeof(schnaps_real));
  schnaps_real  *b = (schnaps_real*) calloc(n,sizeof(schnaps_real));
  schnaps_real  *y = (schnaps_real*) calloc(n,sizeof(schnaps_real));
 
  for (int i=0;i< n; i++)
    xref[i]= ((schnaps_real) rand())/((schnaps_real) RAND_MAX);

  printf("xref: \n");
  for (int i=0;i< n;i++)
    printf("%f\t",xref[i]);
  printf("\n");

  // compute rhs
  MatVect(&lsol,xref,b);

  printf("b: \n");
  for (int i=0;i< n;i++)
    printf("%f\t",b[i]);
  printf("\n");

  lsol.rhs=b;
  lsol.sol=x;

  LUDecompLinearSolver(&lsol);
  DisplayLinearSolver(&lsol);
  SolveLinearSolver(&lsol);

  // compare with reference solution
  for (int i=0;i<n;i++)
    printf("xref:\t %f \t x:\t %f \t diff: \t %f \n",xref[i],x[i],x[i]-xref[i]);

  MatVect(&lsol,x,y);
  printf("b/y: \n");

  for (int i=0;i< n;i++)
    printf("%f\t \t %f\n",b[i],y[i]);

  printf("\n");

  // cleanup
  FreeLinearSolver(&lsol);
}


//
int test_KLU_block_ut(int n, int block_size){
  KLU klumat;
  //
  printf(" KLU test with block triangular matrix\n");
  //
  InitKLU(&klumat, n);
  //
  int *PL = calloc(n,sizeof(int)); // line permutation vector
  int *PC = calloc(n,sizeof(int)); // column permutation vector
  int *PLinv = calloc(n,sizeof(int)); // line permutation vector (inverse)
  int *PCinv = calloc(n,sizeof(int)); // column permutation vector (inverse)

  // set permutation to identity
  for (int i=0;i< n;i++){
    PL[i]=i;
    PC[i]=i;
  }
  // generate random perturbation
  for (int i=0;i<n;i++){
    int r=rand()%n;
    int tmp=PL[r];
    PL[r]=PL[i];
    PL[i]=tmp;
  }  
  for (int i=0;i<n;i++){
    int r=rand()%n;
    int tmp=PC[r];
    PC[r]=PC[i];
    PC[i]=tmp;
  }
  for (int i=0;i< n;i++){
    PLinv[PL[i]]=i;
    PCinv[PC[i]]=i;
  }
  //
  for (int i=0;i < n; i++){
    for (int j=0; j< n; j++){
      if (i/block_size <= j/block_size){
        SwitchOnKLU(&klumat, PL[i], PC[j]);
      }
    }
  }
  //
  printf("Building CSR struct\n");
  AllocateKLU(&klumat);
  DisplayKLU(&klumat);
  printf(" KLU found %i blocks , %i offdiagblocks , largest block size %i\n",
    klumat.symbolic->nblocks,klumat.symbolic->nzoff,klumat.symbolic->maxblock);
  // 
  // cleanup
	free(PL);
	free(PC);
	free(PLinv);
	free(PCinv);
	FreeKLU(&klumat);
}
int Test_KLU_Steady(void) {

  bool test = true;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  schnaps_real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;

  model.m=1; 
  model.NumFlux=KLU_Upwind_NumFlux;
  model.InitData = TestSteady_KLU_InitData; 
  model.ImposedData = TestSteady_KLU_ImposedData; 
  model.BoundaryFlux = KLU_Steady_BoundaryFlux; 
  model.Source = TestSteady_KLU_Source; 

  int deg[]={1, 1, 0};
  int raf[]={3, 3, 1};
  
  CheckMacroMesh(&mesh, deg, raf);

  schnaps_real tmax = 0.001;
  
  Simulation simu2;
  EmptySimulation(&simu2);
  InitSimulation(&simu2, &mesh, deg, raf, &model);
  simu2.dt = tmax / 10;
  ThetaTimeScheme(&simu2, tmax, simu2.dt);
  
  schnaps_real dd = L2error(&simu2);

  printf("erreur implicit L2=%.12e\n", dd);

  PlotFields(0,false, &simu2, "f0", "dgvisu_f0.msh");
  PlotFields(1,false, &simu2, "f1", "dgvisu_f1.msh");
  schnaps_real tolerance = _SMALL;
  test = test && (dd < tolerance);

  FreeMacroMesh(&mesh);
  
  return test;
}



void TestSteady_KLU_ImposedData(const schnaps_real *xy, const schnaps_real t, schnaps_real *w) {

  schnaps_real x=xy[0];
  schnaps_real y=xy[1];


  w[0] = x + y;
  //w[1] = x - y;


}


void TestSteady_KLU_Source(const schnaps_real *xy, const schnaps_real t, const schnaps_real *w, schnaps_real *S){
  
  schnaps_real x=xy[0];
  schnaps_real y=xy[1];
  schnaps_real z=xy[2];

  S[0] = vit[0] * 1 + vit[1] * 1 + vit[2] * 0;
  //S[1] = 0;

}

void TestSteady_KLU_InitData(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  TestSteady_KLU_ImposedData(x, t, w);
}


void KLU_Steady_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux) {
  schnaps_real wR[3];
  TestSteady_KLU_ImposedData(x , t, wR);
  KLU_Upwind_NumFlux(wL, wR, vnorm, flux);
}

void KLU_Upwind_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux){
  schnaps_real flux_temp=0;

  schnaps_real vn = vit[0] * vnorm[0] + vit[1] * vnorm[1] + vit[2] * vnorm[2];
  schnaps_real vnp = vn > 0.0 ? vn : 0.0;
  schnaps_real vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
  //flux[1] = 0;
  
};
