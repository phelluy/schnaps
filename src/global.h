#ifndef _GLOBAL_H
#define _GLOBAL_H
//! brief global variables and defs

#include <stdbool.h>
#ifdef _WITH_STARPU
#include <starpu.h>

#if (STARPU_MAJOR_VERSION >= 1) && (STARPU_MINOR_VERSION >= 2)
#define STARPU_SUPPORT_ARBITER
#else
#warning StarPU Arbiter is not supported
#endif
void init_global_arbiter();
void register_data_arbiter(starpu_data_handle_t handle);
void destroy_global_arbiter();

#endif //_WITH_STARPU

#ifndef _OPENMP
// activate pthread if openmp is not here
//#define _WITH_PTHREAD
#endif
#ifdef _WITH_OPENCL
extern int nplatform_cl;
extern int ndevice_cl;
char numflux_cl_name[1024]; // FIXME: move to field struct.
char cl_buildoptions[1024];
#endif //_WITH_OPENCL

#define __constant
#define __local
#define __private

#ifndef _DOUBLE_PRECISION
#define schnaps_real float
#else
#define schnaps_real double
#endif

#ifndef _DOUBLE_PRECISION
#define  _VERY_SMALL (1e-7)
#define _SMALL (1e-5)
#else
#define _VERY_SMALL (1e-14)
#define _SMALL (1e-10)
#endif

#pragma start_opencl
typedef struct KineticData{
  int time_order;
  int nb_elem_v;
  int deg_v;
  int mv;
  int index_max_kin;
  int index_rho;
  int index_phi;
  int index_ex;
  int index_ey;
  int index_ez;
  int index_u;
  int index_P;
  int index_T;
  int index_max;
  schnaps_real vmax;
  schnaps_real dv;
  schnaps_real gamma;
  schnaps_real knud;
  bool solve_quasineutrality;
  bool substract_mean_charge;
  // quasi neutrality damping (term in front of (phi -phibar) in
  // QN or damped Laplace equation)
  schnaps_real qn_damping;
} KineticData;
#pragma end_opencl


typedef struct LatticeData{
  int q;
  int d;
  schnaps_real ** q_tab;
  schnaps_real * w_tab;
  
  int temp_const;
  schnaps_real c;
  schnaps_real tau;
  int index_max_q;
  int index_rho;
  int index_ux;
  int index_uy;
  int index_uz;
  int index_temp;
  int index_p;
  int index_max;
  // equilibrium distribution function 
  // in : int i : index of velociy node
  // in : void self : pointer to lattice data, should be cast inside the function
  schnaps_real (*feq) (int i_node,void *self,schnaps_real rho,schnaps_real ux,schnaps_real uy,schnaps_real uz,schnaps_real temp,schnaps_real p);
  void (*collect_diags) (void * simu, schnaps_real *diag_vals);
  //
  schnaps_real diag_2d_period;
  // in the case the fi are splitt across multiple simu object, we store here a global node index selector
  // this is used by the Lattice_OneNodeNumFlux to select the proper velocity on the lattice.
  int current_node_index;
  // 
} LatticeData;

//!\bried forward declaration of container for lattice boltzmann model parameters
typedef struct LBModelDescriptor LBModelDescriptor;
typedef struct LBMSimulation LBMSimulation;
//
typedef struct LatticeBoltzmannSimData{
  int current_node_index; //! shared velocity node selector
  int current_advec_node_index; //! shared velocity node selector for advection only
  LBMSimulation *current_lb_sim; //! shared pointer to current simulation (for compat with rk schemes routines);
  LBModelDescriptor *lb_model; //! array[nb_lattices] pointing to individual lattice models. 
} LatticeBoltzmannSimData;
//

extern KineticData schnaps_kinetic_data;
extern LatticeData schnaps_lattice_data;
extern LatticeBoltzmannSimData schnaps_lbm_simdata;

void InitKineticData(KineticData *kd, int nbelemv, int degv);
void InitLatticeData(LatticeData *ld, int dim, int Q,int temp,schnaps_real sound);
void InitLatticeBoltzmannSimData(LatticeBoltzmannSimData *lsd);

extern bool starpu_is_init;
extern bool starpu_use;
// Use of c kernels when use of starpu
extern bool starpu_c_use;
// Use of opencl kernels when use of starpu
extern bool starpu_ocl_use;
// Use of cuda kernels when use of starpu
extern bool starpu_cuda_use;


// compute and get the charge from opencl computationsd
extern bool schnaps_ocl_getcharge;

// OpenCL program for StarPU
extern bool opencl_program_is_init;
#ifdef _WITH_STARPU
extern struct starpu_opencl_program opencl_program;

#define LOAD_OPENCL_PROGRAM_SPU()               \
  if (!opencl_program_is_init) {                \
    opencl_program_is_init = true;              \
    printf("load OpenCL program...\n");         \
    STARPU_CHECK_RETURN_VALUE(                  \
        starpu_opencl_load_opencl_from_file(    \
            "./schnaps.cl",                     \
            &opencl_program,                    \
            cl_buildoptions),                   \
        "starpu_opencl_load_opencl_from_file"); \
    printf("spu_debug: %s \n",cl_buildoptions);	\
  }

#else

#define LOAD_OPENCL_PROGRAM_SPU()               \
  assert(opencl_program_is_init);

#endif //_WITH_STARPU



#endif // #ifndef _GLOBAL_H
