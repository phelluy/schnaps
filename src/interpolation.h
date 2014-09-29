#ifndef _INTERPOLATION_H
#define _INTERPOLATION_H

void hexa_subcell_lobatto_ref_pg_vol(int* param,int ipg,double* xpg,double* wpg);
void hexa_subcell_lobatto_ref_pg_face(int* param,int ifa,int ipg,double* xpg,double* wpg);
void hexa_subcell_lobatto_psi_ref(int* param, int ib, double* xref, double* psi, double* dpsi);
void hexa_subcell_lobatto_psi_ref(int* param, int ib, double* xref, double* psi, double* dpsi);
void hexa_subcell_lobatto_grad_psi_pg(int* param,int ib,int ipg,REAL* dpsi);

npgf(deg,ifa);
npg(deg);


#endif
