#ifndef _LOBATTO_H
#define _LOBATTO_H


// return weight of glop i with i=0..deg
double wglop(int deg, int i);

// return glop i with i=0..deg
double glop(int deg, int i);

// return the 1d derivative of lagrange polynomial ib at glop ipg
double dlag(int deg, int ib, int ipg);


#endif
