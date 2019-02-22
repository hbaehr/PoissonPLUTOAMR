/* ///////////////////////////////////////////////////////////////////// */
/*!                                                                                                                                                                                                                                                                           
  \file                                                                                                                                                                                                                                                                       
  \brief Compute rhs for self-gravity
                                                                                                                                                                                                                         
  Compute the one-dimensional right hand side for the
  self-gravity in the direction given by ::g_dir.
                                                                                                                                                                                                                                      
  \authors H. Baehr (baehr7@gmail.com)\n
                                                                                                                                                                                                                                      
 \b References
                                                                                                                                                                                                                                                          
 \date   February 21, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void SG_RHS (const Data *d, Data_Arr dU, double *dcoeff,
             double **aflux, double dt, int beg, int end, Grid *grid)
/*!                                                                                                                                                                                                                                                                           
 * \param [in]   d           pointer to PLUTO Data structure 
 * \param [out]  dU          a 4D array containing conservative variables                                                                                                                                                                                                     
 *                           increment                                                                                                                                                                                                                                        
 * \param [out]  dcoeff      1D array of diffusion coefficients                                                                                                                                                                                                               
 * \param [out]  aflux       pointer to 2D array for AMR re-fluxing                                                                                                                                                                                                           
 *                           operations                                                                                                                                                                                                                                       
 * \param [in]   dt          the current time-step                                                                                                                                                                                                                            
 * \param [in]   beg,end     initial and final interface indices                                                                                                                                                                                                              
 * \param [in]   grid        pointer to Grid structure.                                                                                                                                                                                                                       
 *                                                                                                                                                                                                                                                                            
 *********************************************************************** */
{
  int i = g_i;
  int j = g_j;
  int k = g_k;
  int nv;
  double dtdV, dtdx;
  static  double *fA;
  static  Sweep sweep;

  /* --------------------------------------------------------                                                                                                                                                                                                                 
   0. Allocate memory                                                                                                                                                                                                                                                         
   -------------------------------------------------------- */

  if (sweep.vn == NULL) {
    MakeState (&sweep);
    fA = ARRAY_1D(NMAX_POINT, double);
  }

  /* --------------------------------------------------------                                                                                                                                                                                                                 
   1. Compute TC flux                                                                                                                                                                                                                                                         
   -------------------------------------------------------- */
  /*
  if (g_dir == IDIR) {
    ITOT_LOOP (i) NVAR_LOOP(nv) sweep.vn[i][nv] = d->Vc[nv][k][j][i];
  } else if (g_dir == JDIR) {
    JTOT_LOOP (j) NVAR_LOOP(nv) sweep.vn[j][nv] = d->Vc[nv][k][j][i];
  } else if (g_dir == KDIR) {
    KTOT_LOOP (k) NVAR_LOOP(nv) sweep.vn[k][nv] = d->Vc[nv][k][j][i];
  }
  SG_Flux (d->Tc, &sweep, dcoeff, beg-1, end, grid);
  */
}
