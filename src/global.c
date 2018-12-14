#include "global.h"

int nplatform_cl = 0;
int ndevice_cl = 0;

bool starpu_is_init = false;
bool starpu_use = false;
bool starpu_c_use = false;
bool starpu_ocl_use = false;
bool starpu_cuda_use = false;
bool schnaps_ocl_getcharge = false;


#pragma start_opencl
#ifndef _NBELEMV 
#define _NBELEMV 2
#define _DEGV 4
#define _VMAX 6.
#endif
#define _MV  (_NBELEMV * _DEGV + 1)
__constant KineticData  schnaps_kinetic_data={
  .nb_elem_v = _NBELEMV,
  .deg_v = _DEGV,
  .mv = _NBELEMV * _DEGV + 1,
  .index_max_kin = _MV - 1,
  .index_rho = _MV ,
  .index_phi = _MV + 1 ,
  .index_ex = _MV + 2,
  .index_ey = _MV + 3,
  .index_ez = _MV + 4,
  .index_u = _MV + 5,
  .index_P = _MV + 6,
  .index_T = _MV + 7,
  .index_max = _MV + 8,
  .vmax = 6.,
  .dv = 2 * _VMAX / _NBELEMV,
  .gamma = 3.,
  .knud =1.,
  .solve_quasineutrality = false,
  .substract_mean_charge = false,
  .qn_damping = 0,
  .time_order=2
};
#pragma end_opencl



LatticeData  schnaps_lattice_data;
// AcousticData schnaps_acoustic_data;
LatticeBoltzmannSimData schnaps_lbm_simdata;
// OpenCL program for StarPU
bool opencl_program_is_init = false;
#ifdef _WITH_STARPU
struct starpu_opencl_program opencl_program;
#endif //_WITH_STARPU

void InitKineticData(KineticData *kd, int nbelemv, int degv){
  kd->nb_elem_v = nbelemv;
  kd->deg_v = degv;
  kd->mv = nbelemv * degv + 1;
  kd->index_max_kin = kd->mv - 1;
  kd->index_rho = kd->mv ;
  kd->index_phi = kd->mv + 1;    
  kd->index_ex = kd->mv + 2;
  kd->index_ey = kd->mv + 3;
  kd->index_ez = kd->mv + 4;
  kd->index_u = kd->mv + 5;
  kd->index_P = kd->mv + 6;
  kd->index_T = kd->mv + 7;
  kd->index_max = kd->mv + 8;
  kd->vmax = 6.;
  kd->dv = 2 * kd->vmax / nbelemv;
  kd->gamma = 3.;
  kd->knud =1.;
  kd->solve_quasineutrality = false;
  kd->substract_mean_charge = false;
  kd->qn_damping = 0;
  kd->time_order=2;
}



void InitLatticeData(LatticeData *ld, int dim, int Q,int temp,schnaps_real sound){

  ld->q=Q;
  ld->d=dim;
  ld->temp_const=temp;
  ld->c=sound;
  ld->index_max_q=Q-1;
  ld->index_rho=Q;
  ld->index_ux=Q+1;
  ld->index_uy=Q+2;
  ld->index_uz=Q+3;
  ld->index_temp=Q+4;
  ld->index_p=Q+5;
  ld->index_max=Q+6;
  //
  ld->q_tab =(schnaps_real **) calloc(Q,sizeof(schnaps_real*));
  ld->w_tab =(schnaps_real *) calloc(Q,sizeof(schnaps_real));
  for(int i=0;i<Q;i++){
    ld->q_tab[i] =(schnaps_real *) calloc(dim,sizeof(schnaps_real));
  }
  // basic D2Q9 lattice nodes
  if ((ld->d== 2) && (ld->q==9)){
  ld->q_tab[0][0]=0.0;
  ld->q_tab[0][1]=0.0;
  ld->q_tab[1][0]=1.0;
  ld->q_tab[1][1]=0.0;
  ld->q_tab[2][0]=0.0;
  ld->q_tab[2][1]=1.0;
  ld->q_tab[3][0]=-1.0;
  ld->q_tab[3][1]=0.0;
  ld->q_tab[4][0]=0.0;
  ld->q_tab[4][1]=-1.0;
  ld->q_tab[5][0]=1.0;
  ld->q_tab[5][1]=1.0;
  ld->q_tab[6][0]=-1.0;
  ld->q_tab[6][1]=1.0;
  ld->q_tab[7][0]=-1.0;
  ld->q_tab[7][1]=-1.0;
  ld->q_tab[8][0]=1.0;
  ld->q_tab[8][1]=-1.0;
  //
  ld->w_tab[0]=4.0/9.0;
  ld->w_tab[1]=1.0/9.0;
  ld->w_tab[2]=1.0/9.0;
  ld->w_tab[3]=1.0/9.0;
  ld->w_tab[4]=1.0/9.0;
  ld->w_tab[5]=1.0/36.0;
  ld->w_tab[6]=1.0/36.0;
  ld->w_tab[7]=1.0/36.0;
  ld->w_tab[8]=1.0/36.0;
  //
  // scale all nodes velocities 
  for (int i=0; i<= ld->index_max_q; i++){
    for (int j=0; j< ld->d;j++){
    schnaps_real utmp=ld->q_tab[i][j]; 
    ld->q_tab[i][j] = utmp * ld->c;
  };
  };
  //
  ld->diag_2d_period=0.0;
  //
  }
}

//
void InitLatticeBoltzmannSimData(LatticeBoltzmannSimData *lsd){
  lsd->current_node_index=0;
  lsd->lb_model=NULL;
  
};

#if defined(_WITH_STARPU)

#if defined(STARPU_SUPPORT_ARBITER)
static starpu_arbiter_t arbiterGlobal;
#endif // defined(STARPU_SUPPORT_ARBITER)
static int arbiterGlobalIsInit = 0;

void init_global_arbiter(){
    assert(arbiterGlobalIsInit == 0);
#if defined(STARPU_SUPPORT_ARBITER)
    arbiterGlobal = starpu_arbiter_create();
#endif // defined(STARPU_SUPPORT_ARBITER)
    arbiterGlobalIsInit = 1;
}
void register_data_arbiter(starpu_data_handle_t handle){
    assert(arbiterGlobalIsInit != 0);
#if defined(STARPU_SUPPORT_ARBITER)
    starpu_data_assign_arbiter(handle, arbiterGlobal);
#endif // defined(STARPU_SUPPORT_ARBITER)
}
void destroy_global_arbiter(){
    assert(arbiterGlobalIsInit != 0);
#if defined(STARPU_SUPPORT_ARBITER)
    starpu_arbiter_destroy(arbiterGlobal);
#endif // defined(STARPU_SUPPORT_ARBITER)
    arbiterGlobalIsInit = 0;
}
#endif // defined(_WITH_STARPU)

