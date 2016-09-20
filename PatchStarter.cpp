#include <cstdio>
#include <string>
using std::string    

/* ***************************************************************************
*  Template for making an FArrayBox to hold the gravitational potential data
*
*  Borrowed from the existing PLUTO file /Src/Chombo/PatchStartup.cpp
*
*  What makes this a little easier than the what is meant for 
*  is the fact that I only need to store one variable, not several
*************************************************************************** */
    
#include "PatchPluto.H"
    
void PatchPluto::starter(FArrayBox& a_gravpot)
{
  int i, j, k;
  int isub, jsub, ksub, nsub = 5;
  int nxtot, nytot, nztot;
  int nv,  l_convert;             // cannot find any purpose for l_convert
  int ibg, ieg, jbg, jeg, kbg, keg;
  //static double **ucons, **uprim;
  double x1,  x2,  x3;
  double x1s, x2s, x3s;
  double dx1, dx2, dx3;
  double gps[256], gp_av[256], b[3];
  double scrh, dx, cylr;

  double ***GP[1];
  Grid   grid[3];
  Box curBox = a_gravpot.box();
  RBox tbox;

  FArrayBox dV(curBox,CHOMBO_NDV);

  curBox.grow(-m_numGhost);

  setGrid(curBox,grid,dV);  /* -- define grid -- */

  jbg = jeg = kbg = keg = 0;

  D_EXPAND(ibg = a_gravpot.loVect()[0]; ieg = a_gravpot.hiVect()[0];  ,
           jbg = a_gravpot.loVect()[1]; jeg = a_gravpot.hiVect()[1];  ,
           kbg = a_gravpot.loVect()[2]; keg = a_gravpot.hiVect()[2];)

  NX1_TOT = nxtot = ieg - ibg + 1;
  NX2_TOT = nytot = jeg - jbg + 1;
  NX3_TOT = nztot = keg - kbg + 1;

 if (uprim == NULL){
   //uprim = ARRAY_2D(NMAX_POINT, NVAR, double);
   //ucons = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

 
 GP[0] = ArrayMap(nztot,nytot,nxtot,a_gravpot.dataPtr(0));

 /* ----  set labels  ---- */

  EXPAND(MXn = VXn = VX1;  ,
         MXt = VXt = VX2;  ,
         MXb = VXb = VX3;)

  #if PHYSICS == MHD || PHYSICS == RMHD
   EXPAND(BXn = BX1;  ,
          BXt = BX2;  ,
          BXb = BX3;)
  #endif

  tbox.ib = 0; tbox.ie = nxtot-1;
  tbox.jb = 0; tbox.je = nytot-1;
  tbox.kb = 0; tbox.ke = nztot-1;
  
/* --------------------------------------------------------------
                    Assign initial conditions
   -------------------------------------------------------------- */

  BOX_LOOP(&tbox,k,j,i){
    x3 = grid[KDIR].x[k];
    x2 = grid[JDIR].x[j];
    x1 = grid[IDIR].x[i];

    GP[0][k][j][i] = gp_av[0] = 0.0;

/*  ----------------------------------------------------------------
                Compute volume averages 
    ---------------------------------------------------------------- */

    #if INITIAL_SMOOTHING == YES

     for (ksub = 0; ksub < nsub; ksub++){
     for (jsub = 0; jsub < nsub; jsub++){
     for (isub = 0; isub < nsub; isub++){

       x1s = x1 + (double)(1.0 - nsub + 2.0*isub)/(double)(2.0*nsub)*dx1;
       x2s = x2 + (double)(1.0 - nsub + 2.0*jsub)/(double)(2.0*nsub)*dx2;
       x3s = x3 + (double)(1.0 - nsub + 2.0*ksub)/(double)(2.0*nsub)*dx3;

       Init (gps, x1s, x2s, x3s);
       for (nv = 0; nv < 1; nv++) {
         gp_av[0] += gps[0]/(double)(nsub*nsub*nsub);
       }
     }}} // ??????

    #else

     Init (gp_av, x1, x2, x3);
    
    #endif

    UU[0][k][j][i] = gp_av[0];

    #if (PHYSICS == MHD || PHYSICS == RMHD) 
     #if ASSIGN_VECTOR_POTENTIAL == YES
      VectorPotentialDiff(b, i, j, k, grid);
      for (nv = 0; nv < DIMENSIONS; nv++) GP[BX1+nv][k][j][i] = b[nv];
     #endif  /* ASSIGN_VECTOR_POTENTIAL */
    #endif /* PHYSICS == MHD || PHYSICS == RMHD */

  }

/* --------------------------------------------------------------------
     Convert primitive variables to conservative ones
   -------------------------------------------------------------------- */

  for (k = 0; k < nztot ; k++) {
  for (j = 0; j < nytot ; j++) {

    for (i = 0; i < nxtot; i++) {
      for (nv = 0 ; nv < 1 ; nv++) {
        us[nv] = uprim[i][nv] = UU[nv][k][j][i];
      }

  /* -- check if primitive values are physically ok -- */

      if (us[RHO] <= 0.0) {
        print ("! startup: density is negative\n");
        QUIT_PLUTO(1);
      }
      #if EOS != ISOTHERMAL
       if (us[PRS] <= 0.0) {
         print ("! startup: pressure is negative\n");
         QUIT_PLUTO(1);
       }
      #endif
      #if (PHYSICS == RHD) || (PHYSICS == RMHD)
       scrh = EXPAND(us[VX1]*us[VX1], + us[VX2]*us[VX2], + us[VX3]*us[VX3]);
       if (scrh >= 1.0){
         print ("! startup: total velocity exceeds 1\n"); 
         QUIT_PLUTO(1);
       }
      #endif
    }

    PrimToCons (uprim, ucons, 0, nxtot-1);
    for (i = 0; i < nxtot; i++) {
    for (nv = 0; nv < 1; nv++) {
      UU[nv][k][j][i] = ucons[i][nv];
    }}

// There should not be any issue between entropy and grav. potential
#if ENTROPY_SWITCH
    Entropy(uprim, UU[ENTR][k][j], 0, nxtot-1);  /* -- primitive: s -- */
    for (i = 0; i < nxtot; i++) {
      UU[ENTR][k][j][i] *= UU[RHO][k][j][i];   /* -- conservative: s*D -- */
    }
#endif

  }}

/* --------------------------------------------------
     Pass GP*dV/m_dx^3 to the library
     
     The passing of gravpot to the library needs to 
     happen, but I am unsure about the necessity 
     of the volume calculation.
   -------------------------------------------------- */

  #if GEOMETRY != CARTESIAN
   #if CHOMBO_CONS_AM == YES
    #if ROTATING_FRAME == YES
     Box aBox = a_gravpot.box();
     for(BoxIterator bit(aBox); bit.ok(); ++bit) {
       const IntVect& iv = bit();
       a_gravpot(iv,iMPHI) += a_gravpot(iv,RHO)*dV(iv,1)*g_OmegaZ;
       a_gravpot(iv,iMPHI) *= dV(iv,1);
     }
    #else
     a_gravpot.mult(dV,1,iMPHI);
    #endif
   #endif
   a_gravpot.mult(dV,0,0);
  #else
   if (g_stretch_fact != 1.) a_gravpot *= g_stretch_fact;
  #endif
 
  FreeArrayMap(GP[0]);
  FreeGrid(grid);
}