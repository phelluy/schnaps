#include "aderdg.h"
#include <assert.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define _MAX(a,b) ((a) > (b) ? (a) : (b))
#define _MIN(a,b) ((a) < (b) ? (a) : (b))

// Stretching fn
#define _STRETCH(x) 0.9*x*x*x+0.1*x

// Final time
#define _TMAX 0.1

// Activate ADER-DG predictor
#define _ADER



int main() {
  ADERDG adg;
  ADERDG_Init(&adg, -1, 1);

  /* for (int ie = 1; ie <= _NBELEMS_IN; ie++) { */
  /*   ADERDG_Predictor(&adg, ie, 0.05); */
  /* } */

  /* adg.tnow = 0.05; */

  /* Plot(&adg); */

  /* assert(1==2); */

  struct timeval start;
  gettimeofday(&start, NULL);

  ADERDG_Solve(&adg, _TMAX);

  struct timeval end;
  gettimeofday(&end, NULL);
  printf("Time: %lu ms\n",
         (end.tv_sec - start.tv_sec) * 1000 +
         (end.tv_usec - start.tv_usec) / 1000);

  ADERDG_Plot(&adg);

  return 0;
}


double stretching(double x) {
  //const double alpha = 2;
  //const double beta = 2;
  // return x;
  //return
  //    alpha * x +
  //    (3 - 2 * alpha - beta) * x * x +
  //    (alpha - 2 + beta) * x * x * x;
  return _STRETCH(x);
}


void ADERDG_Init(ADERDG* adg, double xmin, double xmax) {
  // Mesh boundaries
  adg->xmin = xmin;
  adg->xmax = xmax;
  assert(xmin < xmax);

  // Interface coordinates
  for (int ifa = 0; ifa < _NBFACES; ifa++) {
    // Apply a streching to get different dx
    const double x = stretching((double)ifa / _NBELEMS_IN);
    adg->face[ifa] = xmin + x * (xmax - xmin);
    //printf("face ifa=%d x=%f\n", ifa, adg->face[ifa]);
  }


  // Get maximal cell size
  double min_dx = adg->face[1] - adg->face[0];
  double max_dx = min_dx;
  for (int ie = 1; ie <= _NBELEMS_IN; ie++) {
    const double dx = adg->face[ie] - adg->face[ie - 1];
    min_dx = _MIN(min_dx, dx);
    max_dx = _MAX(max_dx, dx);
  }
  printf("Min dx: %f\n", min_dx);
  printf("Max dx: %f\n", max_dx);
  assert(min_dx > 0);
  assert(max_dx > 0);


  // Get the maximal time step
  adg->dt_small = _CFL * min_dx;  // TODO: put the velocity


  // Compute the levels
  adg->ncfl = 0;
  for (int ie = 1; ie <= _NBELEMS_IN; ie++) {
    const double dx = adg->face[ie] - adg->face[ie - 1];
    adg->cell_level[ie] = (int)(log(dx / min_dx) / log(2));
    adg->ncfl = _MAX(adg->ncfl, adg->cell_level[ie]);
  }
  
  // Impose the maximal number of levels
  //adg->ncfl = 1;
  int* nb_el = (int*) malloc(sizeof(int) * (adg->ncfl + 1));
  for (int i = 0; i <= adg->ncfl; i++) 
    nb_el[i] = 0;
  for (int ie = 0; ie < _NBELEMS; ie++) {
    //if (adg->cell_level[ie] > adg->ncfl)
      //adg->cell_level[ie] = adg->ncfl;
    nb_el[adg->cell_level[ie]]++;
  }
  printf("Found %d cfl levels\n", adg->ncfl + 1);
  for (int i = 0; i <= adg->ncfl; i++) 
    printf("  Elements of level %d: %d\n", i, nb_el[i]);
  free(nb_el);

  // Revert cell levels:
  // 0: smallest cells
  // ncfl: biggest cells
//  for (int ie = 0; ie < _NBELEMS; ie++)
//    adg->cell_level[ie] = adg->ncfl - adg->cell_level[ie];

  // Convention: boundary cells have the same level than their neighbors
  adg->cell_level[0] = adg->cell_level[1];
  adg->cell_level[_NBELEMS_IN + 1] = adg->cell_level[_NBELEMS_IN];

  // Interface levels
  for (int ifa = 0; ifa < _NBFACES; ifa++) {
    adg->face_level[ifa] = _MIN(adg->cell_level[ifa], adg->cell_level[ifa + 1]);

    // Check the level propagation
    assert(adg->cell_level[ifa] - adg->cell_level[ifa + 1] < 2);
    assert(adg->cell_level[ifa] - adg->cell_level[ifa + 1] > -2);
  }

  // Create index of interfaces ordered according to their level (0 -> ncfl)
  int* pos = adg->face_order;
  for (int level = 0; level <= adg->ncfl; level++)
    for (int ifa = 0; ifa < _NBFACES; ifa++)
      if (adg->face_level[ifa] == level)
        *pos++ = ifa;
  assert(pos - adg->face_order == _NBFACES);


  // Get the maximal time step
  adg->dt = adg->dt_small * (1 << adg->ncfl);
  printf("Min dt: %f\n", adg->dt_small);
  printf("Max dt: %f\n", adg->dt);


  // Init internal field data
  // Boundary elements will get their exact field computed at each step
  adg->tnow = 0;
  for (int ie = 0; ie < _NBELEMS; ie++)
    adg->cell_tnow[ie] = 0;
  for (int ie = 1; ie <= _NBELEMS_IN; ie++) {
    for (int ipg = 0; ipg < _NGLOPS; ipg++) {
      const double h = adg->face[ie] - adg->face[ie - 1];
      const double x = adg->face[ie - 1] + h * glop(_D, ipg);
      ExactSol(x, adg->tnow, adg->wnow[ie][ipg]);

      for(int iv = 0; iv < _M; iv++)
	adg->wnext[ie][ipg][iv] = adg->wnow[ie][ipg][iv];
    }
  }


  // Init the computed prediction indicator
  for (int ie = 0; ie < _NBELEMS; ie++)
    adg->pred_done[ie] = false;
}


void ADERDG_Plot(ADERDG* adg) {
  FILE * gnufile;
  gnufile = fopen("adgplot.dat", "w");

  double l2error = 0;
  double min_h = 10;

  for (int ie = 1; ie <= _NBELEMS_IN; ie++) {
    for (int ipg = 0; ipg < _NGLOPS; ipg++) {
      const double h = adg->face[ie] - adg->face[ie - 1];
      const double x = adg->face[ie - 1] + h * glop(_D, ipg);

      double wexact[_M];
      ExactSol(x, adg->tnow, wexact);

      fprintf(gnufile, "%f ", x);  // => 1

      //fprintf(gnufile, "%f ", adg->wnow[ie][ipg][0]);
      //fprintf(gnufile, "%f ", wexact[0]);
      l2error += (adg->wnow[ie][ipg][0] - wexact[0]) * (adg->wnow[ie][ipg][0] - wexact[0]) * wglop(_D, ipg) * h;
      min_h = (h < min_h) ? h : min_h;


      // Numerical values => 2 & 3
      for (int iv = 0; iv < _M; ++iv)
	fprintf(gnufile, "%f ", adg->wnow[ie][ipg][iv]);

      // Exact values => 4 & 5
      for (int iv = 0; iv < _M; ++iv)
	fprintf(gnufile, "%f ", wexact[iv]);

      // Predictor (debug) => 6 & 7
      for (int iv = 0; iv < _M; ++iv)
	fprintf(gnufile, "%f ", adg->wpred[ie][ipg][iv]);

      // Cell levels => 8
      fprintf(gnufile, "%d ", adg->cell_level[ie]);

      // EOL
      fprintf(gnufile, "\n");
    }
  }

  printf("L2 Error: %f\n", sqrt(l2error));
  printf("log: %f\n", log(sqrt(l2error)));
  printf("min dx: %f\n", min_h);
  printf("log: %f\n", log(min_h));

  fclose(gnufile);
}


void ADERDG_Solve(ADERDG* adg, double tmax) {
  int iter = 0;
#ifdef _ADER
  // TODO finalize loop with normal time steps ?
  const int itermax = tmax / adg->dt + 1;
#else
  const int itermax = tmax / adg->dt_small + 1;
#endif
  const int freq = (1 >= itermax / 10) ? 1 : itermax / 10;

  while (iter < itermax) {
#ifdef _ADER
    ADERDG_BigStep(adg);
    adg->tnow += adg->dt;
#else
    ADERDG_TimeStep(adg);
    adg->tnow += adg->dt_small;
#endif
    iter++;

    if (iter % freq == 0)
      printf("iter=%d/%d t=%f dt=%f dt_small=%f\n", iter, itermax, adg->tnow,adg->dt,adg->dt_small);
  }
}


void ADERDG_Predictor(ADERDG* adg, int ie, double tpred) {
  assert(ie < _NBELEMS);

  // Control that prediction is not already computed
  if (adg->pred_done[ie]) return;


  // Impose the exact values on boundary cells
  if (ie == 0) {
    // iv = _D: last point of left cell
    ExactSol(adg->face[0], tpred, adg->wpred[0][_D]);
    return;
  } else if (ie == _NBELEMS - 1) {
    // iv = 0 : first point of right cell
    ExactSol(adg->face[_NBFACES - 1], tpred, adg->wpred[_NBELEMS - 1][0]);
    return;
  }


  const double h = adg->face[ie] - adg->face[ie - 1];

  // Loop over conservative vars
  for (int iv = 0; iv < _M; iv++) {
    // Load field of every gauss point
    double wnow[_NGLOPS];
    for (int ipg = 0; ipg < _NGLOPS; ipg++)
      wnow[ipg] = adg->wnow[ie][ipg][iv];

    // Get the prediction weight
    const double dt = (tpred - adg->cell_tnow[ie]);
    const double pw = dt / h * velocity[iv];
    // printf("ie=%d tpred=%f tnow=%f dt=%f\n", ie, tpred, adg->cell_tnow[ie], dt);

    // Get the prediction
    double wpred[_NGLOPS];
    Prediction(wnow, pw, wpred);

    // Store prediction
    for (int ipg = 0; ipg < _NGLOPS; ipg++)
      adg->wpred[ie][ipg][iv] = wpred[ipg];
  }

  // Update the computed prediction indicator
  adg->pred_done[ie] = true;
}


void ADERDG_TimeStep(ADERDG* adg) {
  // Init the derivatives to zero
  for (int ie = 1; ie <= _NBELEMS_IN; ie++)
    for (int ipg = 0; ipg < _NGLOPS; ipg++)
      for(int iv = 0; iv < _M; iv++)
	adg->dtw[ie][ipg][iv] = 0;


  // Each step is computed on each cell with the smallest dt
  const double dt = adg->dt_small;
  const double tpred = adg->tnow + dt / 2 ;


  // Predict the values of the field at an intermediate time step
  for (int ie = 0; ie < _NBELEMS; ie++) {
    adg->pred_done[ie] = false;
    ADERDG_Predictor(adg, ie, tpred);
  }


  // Compute the face flux terms
  // iv = _D: last point of left cell
  // iv = 0 : first point of right cell
  for (int ifa = 0; ifa < _NBFACES; ifa++) {
    const int ie = ifa;
    double* wL = adg->wpred[ie][_D];
    double* wR = adg->wpred[ie + 1][0];
    double flux[_M];
    NumFlux(wL, wR, flux);
    for (int iv = 0; iv < _M; iv++) {
      adg->dtw[ie][_D][iv] -= flux[iv];
      adg->dtw[ie + 1][0][iv] += flux[iv];
    }
  }


  // Compute the volume terms
  for (int ie = 1; ie <= _NBELEMS_IN; ie++) {
    const double h = adg->face[ie] - adg->face[ie - 1];
    for (int ipg = 0; ipg < _NGLOPS; ipg++) {
      const double omega = wglop(_D, ipg) * h;

      // Flux
      double flux[_M];
      NumFlux(adg->wpred[ie][ipg], adg->wpred[ie][ipg], flux);

      for (int ib = 0; ib < _D + 1; ib++) {
	// Derivative of basis function ib at glop ipg
	const double dpsi = dlag(_D, ib, ipg) / h;
	for (int iv = 0; iv < _M; iv++)
	  adg->dtw[ie][ib][iv] += omega * dpsi * flux[iv];
      }
    }
  }


  // Divide by the mass matrix
  for (int ie = 1; ie <= _NBELEMS_IN; ie++) {
    const double h = adg->face[ie] - adg->face[ie - 1];
    for (int ipg = 0; ipg < _NGLOPS; ipg++) {
      const double omega = wglop(_D, ipg) * h;
      for (int iv = 0; iv < _M; iv++)
	adg->dtw[ie][ipg][iv] /= omega;
    }
  }


  // Update wnext and copy wnext into wnow for the next time step
  for (int ie = 1; ie <= _NBELEMS_IN; ie++) {
    adg->cell_tnow[ie] += dt;
    for (int ipg = 0; ipg < _NGLOPS; ipg++)
      for (int iv = 0; iv < _M; iv++) {
	adg->wnext[ie][ipg][iv] = adg->wnow[ie][ipg][iv] + dt * adg->dtw[ie][ipg][iv];
	adg->wnow[ie][ipg][iv] = adg->wnext[ie][ipg][iv];
      }
  }


  // Update time of first and last cells
  adg->cell_tnow[0] += dt;
  adg->cell_tnow[_NBELEMS_IN + 1] += dt;
}


void ADERDG_BigStep(ADERDG* adg) {
  // Control tnow
  for (int ie = 0; ie < _NBELEMS; ie++) {
    // printf("ie=%d cell_tnow=%f tnow=%f\n", ie, adg->cell_tnow[ie], adg->tnow);
    assert(fabs(adg->cell_tnow[ie] - adg->tnow) < 1e-8);
  }

  // Init the derivatives to zero
  for (int ie = 1; ie <= _NBELEMS_IN; ie++)
    for (int ipg = 0; ipg < _NGLOPS; ipg++)
      for(int iv = 0; iv < _M; iv++)
        adg->dtw[ie][ipg][iv] = 0;


  // Optimized loop
  //  - only 1 step on cells of level 0
  //  - only 2 steps on cells of level 1
  //  - ...
  //  - only 2**ncfl steps on cells of level ncfl
  for (int step = 1; step <= (1 << adg->ncfl); step++) {
    // Get the maximal level to treat for each step
    int max_level = 0;
#ifdef _ADER
    // ADER-DG: only treat cells for which 2**level divides step
    for (int i = 0; i <= adg->ncfl; i++)
      if (step % (1 << i) == 0) max_level = i;
#else
    // Otherwise: treat all the cells with the time step of the smallest cell
    max_level = adg->ncfl;
#endif


    // Compute the face flux terms of cells to treat
    // Predictor is needed for flux and volume
    // On interfaces where cells have not the same level, for the second step of small cell:
    //  - small cell needs prediction of big cell at (tnow + small dt / 2 * 3) for flux term
    //  - big cell needs its own prediction at (tnow + small dt) for volume term
    // In a step, 2 different predictions have to be computed for big cells next to small cells
    // To accomplish that we loop over faces sorted according to their level
    for (int i = 0; i < _NBFACES; i++) {
      const int ifa = adg->face_order[i];

      if (adg->face_level[ifa] > max_level)
        break;

      //if (adg->face_level[ifa] <= max_level)
      {
        const int ieL = ifa;
        const int ieR = ifa + 1;

        // Interface between two levels
        if (adg->cell_level[ieL] != adg->cell_level[ieR]) {
          const int ie_small = (adg->cell_level[ieL] < adg->cell_level[ieR]) ? ieL : ieR;
          const int ie_big = (adg->cell_level[ieL] > adg->cell_level[ieR]) ? ieL : ieR;

          // Flux contributions are computed at the rhythm of small cell
          // Get prediction time of small cell
          const double tpred = adg->tnow +
                               (adg->dt_small / 2) *
                               (2 * step - pow(2, adg->cell_level[ie_small]));

          // Compute predictions
          // Predictor determines the dt to apply according to the level of the cell
          // Predictor applies exact solution on boundary
          ADERDG_Predictor(adg, ieL, tpred);
          ADERDG_Predictor(adg, ieR, tpred);

          // Big cell needs another prediction for its volume term
          // Predicton will be recomputed thanks to interface ordering
          adg->pred_done[ie_big] = false;

          // Compute flux
          double* wL = adg->wpred[ieL][_D];
          double* wR = adg->wpred[ieR][0];
          double flux[_M];
          NumFlux(wL, wR, flux);
          for (int iv = 0; iv < _M; iv++) {
            adg->dtw[ieL][_D][iv] -= flux[iv] / (1 + (ieL == ie_big));
            adg->dtw[ieR][0][iv] += flux[iv] / (1 + (ieR == ie_big));
          }

        } else {
          // Get prediction time of current cells of same level
          const double tpred = adg->tnow +
                               (adg->dt_small / 2) *
                               (2 * step - pow(2, adg->face_level[ifa]));

          // Compute predictions
          // Predictor determines the dt to apply according to the level of the cell
          // Predictor applies exact solution on boundary
          ADERDG_Predictor(adg, ieL, tpred);
          ADERDG_Predictor(adg, ieR, tpred);

          // Compute flux
          double* wL = adg->wpred[ieL][_D];
          double* wR = adg->wpred[ieR][0];
          double flux[_M];
          NumFlux(wL, wR, flux);
          for (int iv = 0; iv < _M; iv++) {
            adg->dtw[ieL][_D][iv] -= flux[iv];
            adg->dtw[ieR][0][iv] += flux[iv];
          }
        }
      }
    }


    // Compute the volume terms of cells to treat
    for (int ie = 1; ie <= _NBELEMS_IN; ie++)
      if (adg->cell_level[ie] <= max_level) {
        const double h = adg->face[ie] - adg->face[ie - 1];
        for (int ipg = 0; ipg < _NGLOPS; ipg++) {
          const double omega = wglop(_D, ipg) * h;

          // Flux
          double flux[_M];
          NumFlux(adg->wpred[ie][ipg], adg->wpred[ie][ipg], flux);

          for (int ib = 0; ib < _D + 1; ib++) {
            // Derivative of basis function ib at glop ipg
            const double dpsi = dlag(_D, ib, ipg) / h;
            for (int iv = 0; iv < _M; iv++)
              adg->dtw[ie][ib][iv] += omega * dpsi * flux[iv];
          }
        }
      }


    // Divide by the mass matrix for the cells to treat
    for (int ie = 1; ie <= _NBELEMS_IN; ie++)
      if (adg->cell_level[ie] <= max_level) {
        const double h = adg->face[ie] - adg->face[ie - 1];
        for (int ipg = 0; ipg < _NGLOPS; ipg++) {
          const double omega = wglop(_D, ipg) * h;
          for (int iv = 0; iv < _M; iv++)
            adg->dtw[ie][ipg][iv] /= omega;
        }
      }


    // Update wnext and copy wnext into wnow for the next time step
    for (int ie = 0; ie < _NBELEMS; ie++)
      if (adg->cell_level[ie] <= max_level) {
        // Get time step of cell
        const double dt = adg->dt_small * pow(2, adg->cell_level[ie]);

        // Update the tnow of treated cells
        adg->cell_tnow[ie] += dt;

        for (int ipg = 0; ipg < _NGLOPS; ipg++) {
          for (int iv = 0; iv < _M; iv++) {
            adg->wnext[ie][ipg][iv] = adg->wnow[ie][ipg][iv] + dt * adg->dtw[ie][ipg][iv];
            adg->wnow[ie][ipg][iv] = adg->wnext[ie][ipg][iv];

            // Nullify dtw for next step
            adg->dtw[ie][ipg][iv] = 0;

            // Reinit prediction indicator for next step
            adg->pred_done[ie] = false;
          }
        }
      }
  }
}



// Model
// -----------------------------------------------------------------------------

void ExactSol(double x, double t, double w[_M]) {
  const double pi = 4 * atan(1.);

  w[0] = cos(2 * pi * (x - t));
  w[1] = cos(2 * pi * (x + t));

  /* w[0] *= w[0] * w[0]; */
  /* w[1] *= w[1] * w[1]; */
}


void Prediction(double* wnow, double pw, double* wpred) {
  switch(_D) {
    case 0:
      wpred[0] = wnow[0];

      //wpred[1] = wnow[1];

      /// etc...
      
      break;
    case 1:
      wpred[0] = (1 + pw) * wnow[0] -
                 pw * wnow[1];

      wpred[1] = pw * wnow[0] +
                 (1 - pw) * wnow[1];

      break;

    case 2:
      wpred[0] = (1 + 3 * pw + 2 * pw * pw) * wnow[0] +
                 (-4 * pw * pw - 4 * pw) * wnow[1] +
                 (2 * pw * pw + pw) * wnow[2];

      wpred[1] = (2 * pw * pw + pw) * wnow[0] +
                 (1 - 4 * pw * pw) * wnow[1] +
                 (2 * pw * pw - pw) * wnow[2];

      wpred[2] = (2 * pw * pw - pw) * wnow[0] +
                 (-4 * pw * pw + 4 * pw) * wnow[1] +
                 (1 + 2 * pw * pw - 3 * pw) * wnow[2];

      break;

    case 3:
      wpred[0] = 0.1000000000e1 * (0.1e1 + 0.6e1 * pw + 0.10e2 * pw * pw + 0.5e1 * pow(pw, 0.3e1)) * wnow[0] - 0.1809016994e1 * (0.1065247584e2 * pw + 0.4472135954e1 + 0.618033988e1 * pw * pw) * pw * wnow[1] + 0.6909830056e0 * (0.2065247584e2 * pw + 0.4472135954e1 + 0.1618033988e2 * pw * pw) * pw * wnow[2] - 0.1000000000e1 * pw * (0.1e1 + 0.5e1 * pw + 0.5e1 * pw * pw) * wnow[3];

      wpred[1] = 0.5000000000e0 * pw * (0.1170820393e2 * pw + 0.3236067977e1 + 0.10e2 * pw * pw) * wnow[0] - 0.3618033987e0 * (-0.2763932023e1 + 0.2763932023e2 * pw * pw + 0.3090169942e2 * pow(pw, 0.3e1)) * wnow[1] + 0.1000000000e1 * pw * (0.1118033988e2 * pw * pw + 0.5e1 * pw - 0.2236067977e1) * wnow[2] - 0.5000000000e0 * pw * (0.1708203931e1 * pw - 0.1236067977e1 + 0.10e2 * pw * pw) * wnow[3];

      wpred[2] = -0.5000000000e0 * pw * (0.1708203931e1 * pw + 0.1236067977e1 - 0.10e2 * pw * pw) * wnow[0] - 0.1000000000e1 * pw * (0.1118033988e2 * pw * pw - 0.5e1 * pw - 0.2236067977e1) * wnow[1] + 0.1381966011e0 * (0.7236067977e1 - 0.7236067977e2 * pw * pw + 0.8090169942e2 * pow(pw, 0.3e1)) * wnow[2] + 0.5000000000e0 * pw * (0.1170820393e2 * pw - 0.3236067977e1 - 0.10e2 * pw * pw) * wnow[3];

      wpred[3] = 0.1000000000e1 * pw * (0.1e1 - 0.5e1 * pw + 0.5e1 * pw * pw) * wnow[0] - 0.1809016993e1 * pw * (0.1708203931e1 + 0.618033988e1 * pw * pw - 0.788854382e1 * pw) * wnow[1] + 0.6909830058e0 * pw * (0.1170820393e2 + 0.1618033988e2 * pw * pw - 0.2788854382e2 * pw) * wnow[2] - 0.1000000000e1 * (-0.1e1 + 0.6e1 * pw - 0.10e2 * pw * pw + 0.5e1 * pow(pw, 0.3e1)) * wnow[3];

      break;

    case 4:
      wpred[0] = 0.9999999983e0 * (0.1e1 + 0.10e2 * pw + 0.14e2 * pow(pw, 0.4e1) + 0.35e2 * pow(pw, 0.3e1) + 0.30e2 * pw * pw) * wnow[0] + 0.6756502479e1 * (-0.2e1 - 0.8417424305e1 * pw - 0.1125227292e2 * pw * pw - 0.4834848610e1 * pow(pw, 0.3e1)) * pw * wnow[1] + 0.5333333325e1 * pw * (0.1e1 + 0.8e1 * pw + 0.7e1 * pow(pw, 0.3e1) + 0.14e2 * pw * pw) * wnow[2] - 0.1410164176e1 * (0.2e1 + 0.1758257570e2 * pw + 0.3874772708e2 * pw * pw + 0.2316515139e2 * pow(pw, 0.3e1)) * pw * wnow[3] + 0.9999999983e0 * pw * (0.1e1 + 0.9e1 * pw + 0.14e2 * pow(pw, 0.3e1) + 0.21e2 * pw * pw) * wnow[4];

      wpred[1] = 0.7142857131e-1 * pw * (0.3474772708e2 + 0.2012340896e3 * pw + 0.3546242389e3 * pw * pw + 0.196e3 * pow(pw, 0.3e1)) * wnow[0] + 0.1378878057e0 * (0.725227292e1 - 0.1692197014e3 * pw * pw - 0.387731045e3 * pow(pw, 0.3e1) - 0.2369075819e3 * pow(pw, 0.4e1)) * wnow[1] + 0.7619047607e0 * pw * (0.14e2 * pw + 0.49e2 * pow(pw, 0.3e1) - 0.4582575695e1 + 0.6415605973e2 * pw * pw) * wnow[2] - 0.2014520251e0 * (-0.7582575695e1 + 0.1158257570e2 * pw + 0.1592340896e3 * pw * pw + 0.1621560597e3 * pow(pw, 0.3e1)) * pw * wnow[3] + 0.7142857131e-1 * pw * (-0.725227292e1 + 0.876591040e1 * pw + 0.1586242389e3 * pw * pw + 0.196e3 * pow(pw, 0.3e1)) * wnow[4];

      wpred[2] = 0.2499999996e0 * pw * (-0.3e1 - 0.6e1 * pw + 0.56e2 * pow(pw, 0.3e1) + 0.28e2 * pw * pw) * wnow[0] + 0.1689125620e1 * (0.1582575695e1 + 0.4834848610e1 * pw - 0.633030278e1 * pw * pw - 0.1933939444e2 * pow(pw, 0.3e1)) * pw * wnow[1] + 0.3333333329e0 * (0.3e1 + 0.112e3 * pow(pw, 0.4e1) - 0.40e2 * pw * pw) * wnow[2] - 0.3525410440e0 * (0.7582575695e1 - 0.2316515139e2 * pw - 0.3033030278e2 * pw * pw + 0.9266060556e2 * pow(pw, 0.3e1)) * pw * wnow[3] + 0.2499999996e0 * pw * (0.3e1 - 0.6e1 * pw + 0.56e2 * pow(pw, 0.3e1) - 0.28e2 * pw * pw) * wnow[4];

      wpred[3] = -0.7142857131e-1 * pw * (-0.725227292e1 - 0.876591040e1 * pw + 0.1586242389e3 * pw * pw - 0.196e3 * pow(pw, 0.3e1)) * wnow[0] + 0.9652146398e0 * (-0.1582575695e1 - 0.2417424305e1 * pw + 0.3323408960e2 * pw * pw - 0.3384394027e2 * pow(pw, 0.3e1)) * pw * wnow[1] - 0.7619047607e0 * pw * (-0.14e2 * pw - 0.49e2 * pow(pw, 0.3e1) - 0.4582575695e1 + 0.6415605973e2 * pw * pw) * wnow[2] - 0.2877886073e-1 * (-0.3474772708e2 + 0.8107802986e3 * pw * pw - 0.1857731045e4 * pow(pw, 0.3e1) + 0.1135092418e4 * pow(pw, 0.4e1)) * wnow[3] - 0.7142857131e-1 * pw * (0.3474772708e2 - 0.2012340896e3 * pw + 0.3546242389e3 * pw * pw - 0.196e3 * pow(pw, 0.3e1)) * wnow[4];

      wpred[4] = 0.9999999983e0 * pw * (-0.1e1 - 0.21e2 * pw * pw + 0.9e1 * pw + 0.14e2 * pow(pw, 0.3e1)) * wnow[0] - 0.1166666665e1 * pw * (-0.2417424305e1 + 0.2125227292e2 * pw - 0.4683484861e2 * pw * pw + 0.28e2 * pow(pw, 0.3e1)) * wnow[1] + 0.5333333325e1 * pw * (-0.1e1 + 0.8e1 * pw + 0.7e1 * pow(pw, 0.3e1) - 0.14e2 * pw * pw) * wnow[2] + 0.1166666665e1 * pw * (0.1158257570e2 - 0.4874772708e2 * pw + 0.6516515139e2 * pw * pw - 0.28e2 * pow(pw, 0.3e1)) * wnow[3] + 0.9999999983e0 * (0.1e1 - 0.10e2 * pw + 0.30e2 * pw * pw - 0.35e2 * pow(pw, 0.3e1) + 0.14e2 * pow(pw, 0.4e1)) * wnow[4];

      break;

    default:
      assert(false);
  }
}


void NumFlux(double* wL, double* wR, double* flux) {
  double vn = velocity[0];
  double vnp = vn > 0 ? vn : 0;
  double vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];

  vn  = velocity[1];
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[1] = vnp * wL[1] + vnm * wR[1];

  /* flux[0] = vn * (wL[0] + wR[0]) / 2;  */
  /* flux[1] = vn * (wL[1] + wR[1]) / 2;  */
}



// Interpolation
// -----------------------------------------------------------------------------

// Gauss LObatto Points (GLOP) up to order 4
static const double gauss_lob_point[] = {
  0.5,
  0,
  1,
  0,
  0.5,
  1,
  0,
  0.276393202250021030359082633127,
  0.723606797749978969640917366873,
  1,
  0,
  0.172673164646011428100853771877,
  0.5,
  0.827326835353988571899146228123,
  1
};

// GLOP weights up to order 4
static const double gauss_lob_weight[] = {
  1,
  0.5,
  0.5,
  0.166666666666666666666666666667,
  0.666666666666666666666666666668,
  0.166666666666666666666666666667,
  0.0833333333333333333333333333333,
  0.416666666666666666666666666666,
  0.416666666666666666666666666666,
  0.0833333333333333333333333333333,
  0.05,
  0.272222222222222222222222222223,
  0.355555555555555555555555555556,
  0.272222222222222222222222222219,
  0.05
};

// Offset for finding the GLOP data for a given degree in the previous arrays
static const int gauss_lob_offset[] = {0, 1, 3, 6, 10};

// GLOP gradient of basis function up to order 4
static const double gauss_lob_dpsi[] = {
  0,
  -1,
  -1,
  1,
  1,
  -3,
  -1,
  1,
  4,
  0,
  -4,
  -1,
  1,
  3,
  -6,
  -1.61803398874989484820458683436,
  .618033988749894848204586834362,
  -1,
  8.09016994374947424102293417177,
  0,
  -2.23606797749978969640917366872,
  3.09016994374947424102293417184,
  -3.09016994374947424102293417182,
  2.23606797749978969640917366872,
  0,
  -8.09016994374947424102293417177,
  1,
  -.618033988749894848204586834362,
  1.61803398874989484820458683436,
  6,
  -10,
  -2.48198050606196571569743868436,
  .75,
  -.518019493938034284302561315632,
  1,
  13.5130049774484800076860550594,
  0,
  -2.67316915539090667050969419631,
  1.52752523165194666886268239794,
  -2.82032835588485332564727827404,
  -5.33333333333333333333333333336,
  3.49148624377587810025755976667,
  0,
  -3.49148624377587810025755976662,
  5.33333333333333333333333333336,
  2.82032835588485332564727827399,
  -1.52752523165194666886268239791,
  2.67316915539090667050969419635,
  0,
  -13.5130049774484800076860550594,
  -1,
  .518019493938034284302561315631,
  -.75,
  2.48198050606196571569743868437,
  10
};

// Offset for finding the GLOP data for a given degree in the previous arrays
static const int gauss_lob_dpsi_offset[] = {0, 1, 5, 14, 30};

double glop(int deg, int ipg) {
  return gauss_lob_point[gauss_lob_offset[deg] + ipg];
}

double wglop(int deg, int ipg) {
  return gauss_lob_weight[gauss_lob_offset[deg] + ipg];
}

double dlag(int deg,int ib,int ipg) {
  return gauss_lob_dpsi[gauss_lob_dpsi_offset[deg] + ib * (deg + 1) + ipg];
}
