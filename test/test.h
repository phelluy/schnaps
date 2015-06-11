#ifndef _TEST_H
#define _TEST_H


// some simple tests
int TestMacroMesh(void);

int TestGeometry(void);

int TestInterpolation(void);

int TestModel(void);

int TestField(void);

int TestCache(void);

int TestFieldDG(void);

int TestFieldRK2(void);

int Test2DMeshDetection(void);

int TestFieldRK2_2D(void);

int TestFieldSubCellDGVol(void);

int TestFieldRK2_2D_SubCell(void);

int TestmEq2(void);

int TestMaxwell2D(void);

int TestGyro(void);

int TestMHD(int argc, char *argv[]);

int TestMHD1D(int argc, char *argv[]);

int TestPIC(void);

int TestPICAccumulate(void);

int TestPeriodic(void);

int TestSkyline(void);

int TestLinearSolver(void);

int TestPoisson(void);

int TestPoisson2d(void);

int Test_TransportVP(void);

int TestCLInfo(void);

int TestKernel(void);

int TestKernelVolume(void);

int TestKernelFlux(void);

int TestKernelInterface(void);

int TestFieldRK2_CL(void);

int TestLandau_Damping_1D(void);

#endif
