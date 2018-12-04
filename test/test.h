#ifndef _TEST_H
#define _TEST_H

#include "schnaps.h"
#include <time.h>


#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

//useless header !!!!!!!!!!!!!!!!!!!
// some simple tests
/* int TestMacroMesh(void); */

/* int TestIGraph(void); */

/* int TestGeometry(void); */

/* int TestInterpolation(void); */

/* int TestModel(void); */

/* int TestField(void); */

/* int TestSimulation(void); */

/* int TestCache(void); */

/* int TestFieldDG(void); */

/* int TestFieldRK2(void); */

/* int Test2DMeshDetection(void); */

/* int TestPCWave(void); */

/* int TestFieldRK2_2D(void); */

/* int TestFieldSubCellDGVol(void); */

/* int TestFieldRK2_2D_SubCell(void); */

/* int TestmEq2(void); */

/* int TestMaxwell2D(void); */

/* int TestGyro(void); */

/* int TestMHD(int argc, char *argv[]); */

/* int TestMHD1D(int argc, char *argv[]); */

/* int TestPIC(void); */

/* int TestPICAccumulate(void); */

/* int TestPeriodic(void); */

/* int TestSkyline(void); */

/* int TestSkyline_SPU(void); */

/* int TestUmfPack(void); */

/* int TestLinearSolver(void); */

/* int TestNonLinearSolver(void); */

/* int TestPoisson(void); */

/* int TestPoisson2d(void); */

/* int Test_TransportVP(void); */

/* int TestCLInfo(void); */

/* int TestKernel(void); */

/* int TestKernelVolume(void); */

/* int TestKernelFlux(void); */

/* int TestKernelInterface(void); */

/* int TestFieldRK2_CL(void); */

/* int TestLandau_Damping_1D(void); */

/* int TestLandauCollision_1D(void); */

/* int Test_Wave_Periodic(void); */

/* int Test_Wave_Steady(void); */

/* int Test_Transport_Steady(void); */

/* int Test_Local_Implicit(void); */

/* int Test_Local_Implicit_SPU(void); */

/* int Test_Graph_Implicit_SPU(void); */

/* int Test_NoMatrix_Implicit_SPU(void); */

/* int Test_Transport_ExImp(void); */

/* int Test_SH_equilibrium(void); */

/* int Test_SH_periodic(void); */

/* int Test_SH_SteadyState_UImposed(void); */

/* int Test_SH_SteadyState_PImposed(void); */

/* int Test_SH_equilibrium_Implicit(void); */

/* int TestOrszagTang(int argc, char *argv[]); */

/* int TestReconnexion(int argc, char *argv[]); */

/* int TestKelvinHelmotz(int argc, char *argv[]); */

/* int TestDoubleTearing(int argc, char *argv[]); */

/* int Test_Wave_Periodic(void); */

/* int Testrealpc(void); */

/* int Test_OrderAdaptivity(void); */

// schnaps_real seconds() {
//   struct timespec ts;
//   schnaps_real res;
// #ifdef __MACH__
//   clock_serv_t cclock;
//   mach_timespec_t mts;
//   host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
//   clock_get_time(cclock, &mts);
//   mach_port_deallocate(mach_task_self(), cclock);
//   ts.tv_sec = mts.tv_sec;
//   ts.tv_nsec = mts.tv_nsec;
//   res =  (schnaps_real)ts.tv_sec + 1e-9 * (schnaps_real)ts.tv_nsec;
// #else
//   clock_gettime(CLOCK_MONOTONIC, &ts);
//   res = (schnaps_real)ts.tv_sec + 1e-9 * (schnaps_real)ts.tv_nsec;
// #endif
//   return res;
// }


#endif
