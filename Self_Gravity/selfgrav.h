/* //////////////////////////////////////////////////////////// */
/* !
   \file
   \brief Self-gravity (SELFGRAV) module header file.
   
   Contains prototypes for the self-gravity module.

   \authors H. Baehr (baehr7@gmail.com)\n

   \date    February 20, 2019
*/
/* //////////////////////////////////////////////////////////// */

void SG_Flux (double ***, const sweep *, double *, int, int, Grid *);
void SG_RHS (const Data *, Data_Arr, double *, double **, double, int, int, Grid *);
void GetGradient(double *** double, **, int, int, Grid *);
