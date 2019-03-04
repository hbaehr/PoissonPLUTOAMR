/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute the force due to self-gravity.

  Compute the acceleration due to to self-gravity along one row of 
  computational zones for the HD and MHD  modules

  \authors H. Baehr (baehr7@gmail.com)\n

  \date    March 4, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void SG_Flux (double ***Phi, const Sweep *sweep,  
              double *dcoeff, int beg, int end, Grid *grid)
/*! 
 * Compute the gradient of the gravitational potential Phi.
 *
 * \param [in]     Phi     3D array containing the dimensionless 
 *                         potential
 * \param [in,out] sweep   pointer to a Sweep structure
 * \param [out]    dcoeff  the diffusion coefficient needed for computing   !!!
 *                         the time step.
 * \param [in]     beg     initial index of computation
 * \param [in]     end     final   index of computation
 * \param [in]     grid    pointer to an array of Grid structures
 *
 * \return This function has no return value.                       
 *********************************************************************** */
{
  int  i, j, k, nv;
  double bgradPhi, dPhi_mag;
  double Fc, Fsat, g1 = g_gamma - 1.0;
  double x1, x2, x3;
  double dvp, dvm, dvl;
  double alpha, uL, uR, bn;
  double vi[NVAR], phi;
  double bck_fld[3];
  double **vc = sweep->vn;
  static double **gradPhi, *pp, *pm;

/* -----------------------------------------------------------
   1. Allocate memory, compute temperature gradient in the
      required direction and set 2 out of the 3 coordinates.
   ----------------------------------------------------------- */

  if (gradT == NULL) {
    gradT = ARRAY_2D(NMAX_POINT, 3, double);
    pp    = ARRAY_1D(NMAX_POINT, double);
    pm    = ARRAY_1D(NMAX_POINT, double);
  }

  GetGradient (T, gradT, beg, end, grid);
  if (g_dir == JDIR || g_dir == KDIR) x1 = grid->x[IDIR][g_i];
  if (g_dir == IDIR || g_dir == KDIR) x2 = grid->x[JDIR][g_j];
  if (g_dir == IDIR || g_dir == JDIR) x3 = grid->x[KDIR][g_k];

/* ----------------------------------------------- 
   2. Compute Thermal Conduction Flux (tcflx).
   ----------------------------------------------- */

  for (i = beg; i <= end; i++){
    
  /* -- 3a. Compute interface values -- */

    NVAR_LOOP(nv) vi[nv]  = 0.5*(vc[i][nv] + vc[i+1][nv]);    

  /* -- 3b. Get TC coefficients at cell interfaces -- */

    if (g_dir == IDIR) x1 = grid->xr[IDIR][i];
    if (g_dir == JDIR) x2 = grid->xr[JDIR][i];
    if (g_dir == KDIR) x3 = grid->xr[KDIR][i];

    dPhi_mag  = D_EXPAND(  gradPhi[i][0]*gradPhi[i][0], 
                       + gradPhi[i][1]*gradPhi[i][1], 
                       + gradPhi[i][2]*gradPhi[i][2]);
    dPhi_mag = sqrt(dPhi_mag) + 1.e-12;

    Acc = dPhi_mag * RHO;

    sweep->sg_flux[i][ENG] = Acc;
  }
}
