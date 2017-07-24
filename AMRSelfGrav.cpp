#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::ifstream;
using std::ios;

// #define HALEM_PROC_SPEED
#ifdef HALEM_PROC_SPEED
#include <cstdio>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#endif

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "FABView.H"

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "memusage.H"

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BCFunc.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"

#include "AMR.H"
#include "AMRLevelPlutoFactory.H"
#include "AMRPoissonOp.H"
#include "AMRMultiGrid.H"

#include "UsingNamespace.H"

#ifdef CH_LINUX
// Should be undefined by default
#define TRAP_FPE
#undef  TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

  // initialize some temporary storage for the gravitational potential
  {
  //  CH_TIME("setup::Udefine"); // This is used to time processes for debugging purposes
    m_gravpot.define(m_grids,m_numCons,m_numGhost*IntVect::Unit);
  }

// Other things which need to be addressed:
// 1) Boundary conditions
// 2) Initial conditions
// 2) Fill the ghost zones of m_gravpot or a_gravpot by interpolation? see LevelPluto for more info

/* Set the right hand side of Poisson equation
*  \nabla \Phi = 4*\pi*G*\rho
*
*  Since everything on the RHS outside of \rho is constant, these are put into
*  the constant value b. However, since \rho has already been calculate in the
*  previous step, it might be unnecessary to do this step.
*
*  I will include it for now for completeness and perhaps testing purposes
*/
void setRHS(Vector<LevelData<FArrayBox>* > a_rhs, // Output array of \rho
            Vector<ProblemDomain>& a_domain,      // Grid domain
            Vector<int>& m_ref_ratio,             // Refinement ratios between levels
            Vector<Real>& a_amrDx,                // *** dx: not sure what this value is atm
            int a_finestLevel)                    // *** number of most refined level
{
  CH_TIME("setRHS");

  for (int lev=0; lev<=a_finestLevel; lev++)      // Looping over all levels
    {
      LevelData<FArrayBox>& levelRhs = *(a_rhs[lev]);                      // Each a_rhs[lev] is an FArrayBox
      const DisjointBoxLayout& levelGrids = levelRhs.getBoxes();           //

      // rhs is cell-centered...
      RealVect ccOffset = 0.5*a_amrDx[lev]*RealVect::Unit;                 //

      DataIterator levelDit = levelGrids.dataIterator();                   // Write data to a level
      for (levelDit.begin(); levelDit.ok(); ++levelDit)                    // Iterate over all points on the level?
        {
          FArrayBox& thisRhs = levelRhs[levelDit];                         // Dummy thisRhs for the iteration of this loop

          if (s_probtype == zeroRHS)                                       // Various patterns for the distribution of RHS
            {
              thisRhs.setVal(0.0);                                         // 0 everywhere
            }
          else if (s_probtype == unityRHS)
            {
              thisRhs.setVal(1.0);                                         // 1 everywhere
            }
          else if (s_probtype == sinusoidal)
            {

              BoxIterator bit(thisRhs.box());                              // Sine wave
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  RealVect loc(iv);
                  loc *= a_amrDx[lev];
                  loc += ccOffset;
                  loc *= 2.0*Pi;

                  thisRhs(iv,0) = D_TERM(sin(loc[0]),
                                         *cos(loc[1]),
                                         *sin(loc[2]));

                 }

            }
          else if (s_probtype == gaussians)
            {
              int numGaussians = 3;                                        // Gaussian distributions
              Vector<RealVect> center(numGaussians,RealVect::Zero);
              Vector<Real> scale(numGaussians, 1.0);
              Vector<Real> strength(numGaussians, 1.0);

              for (int n=0; n<numGaussians; n++)
                {
                  if (n==0)
                    {
                      strength[0] = 1.0;
                      scale[0] = 1.0e-2;
                      center[0] = 0.25*RealVect::Unit;
                    }
                  else if (n == 1)
                    {
                      strength[1] = 3.0;
                      scale[1] = 1.0e-2;
                      center[1] = RealVect(D_DECL(0.5,0.75, 0.75));
                    }
                  else if (n == 2)
                    {
                      strength[2] = 2.0;
                      scale[2] = 1.0e-2;
                      center[2] = RealVect(D_DECL(0.75,0.5, 0.5));
                    }
                  else
                    {
                      MayDay::Error("too many Gaussian sources attempted");
                    }
                }

              thisRhs.setVal(0.0);

              BoxIterator bit(thisRhs.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  RealVect loc(iv);
                  loc *= a_amrDx[lev];
                  loc += ccOffset;

                  for (int n=0; n<numGaussians; n++)
                    {
                      RealVect dist = loc - center[n];
                      Real radSqr = D_TERM(dist[0]*dist[0],
                                           +dist[1]*dist[1],
                                            +dist[2]*dist[2]);

                       Real val = strength[n]*exp(-radSqr/scale[n]);
                       thisRhs(iv,0) += val;
                     }
                 }
             }
           else
             {
               MayDay::Error("undefined problem type");
             }
         } // end loop over grids on this level
     } // end loop over levels
 }







/* ********************************************************** */
void solveSelfGravPot(Vector<LevelData<FArrayBox>* >& a_gravpot,       // Output self-gravity potential: m_gravpot
                      const Vector<LevelData<FArrayBox>* > a_U,  // Input density: m_UNew
                      const Vector<DisjointBoxLayout>& a_grids,    // Grid geometries at all levels:
                      const Vector<int>& m_ref_ratios,              // Refinement ratios between levels: m_ref_ratio
                      const ProblemDomain& a_domain,         // Coarsest domain: m_domain
                      Real alpha, Real beta, Real a_dx)    // constants alpha=0, beta=1/(4*pi*G)
/*
 * Example from the Chombo documentation:
 *
 * solveSelfGravPot solves (alpha I + beta Laplacian) phi = rhs
 * using AMRMultiGrid and AMRPoissonOp
 * Inputs:
 *  rhs: Right-hand side of the solve over the level.
 *  grids: AMRHierarchy of grids
 *  refRatio: refinement ratios
 *  level0Domain: domain at the coarsest AMR level
 *  coarsestDx: grid spacing at the coarsest level
 *  alpha: identity coefficient
 *  beta: Laplacian coefficient
 * Outputs:
 *  phi = (alpha I + beta Lapl)^{-1}(rhs)
 *
 ************************************************************ */
{
int numlevels = a_U.size(); // A different array for each refinement level

//define the operator factory
AMRPoissonOpFactory opFactory;
opFactory.define(m_domain,
                 m_grids, m_ref_ratios, m_dx,
                 &ParseBC, alpha, beta);

//this is the solver we shall use
AMRMultiGrid<LevelData<FArrayBox> > solver;

//this is the solver for the bottom of the muligrid v-cycle (???)
BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;

//bicgstab can be noisy
bottomSolver.m_verbosity = 0;

//define the solver
solver.define(m_domain, opFactory, &bottomSolver, numlevels);

//we want to solve over the whole hierarchy
int lbase = 0;

//so solve already.
solver.solve(m_gravpot, m_U, numlevels-1, lbase);
}

// Loop over the patches of a level to assign initial conditions
void PatchPluto::initiate(LevelData<FArrayBox>& a_gravpot)
{

 CH_assert(m_isDefined);

// DataIterator does what?
 DataIterator dit = a_gravpot.boxLayout().dataIterator();

 // Iterator for all grids in this level
 for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_gravpot[dit()];   // Not U, but what?
    // Set up initial conditions in this patch
    starter(U);                        // NOt U, but what?
  }
}

/* Calling procedure should look like:
 void solveSelfGravPot(m_gravpot,            // Output self-gravity potential
                       m_U,                  // Input density: m_UNew
                       m_grids,              // Grid geometries
                       m_ref_ratio,          // Vector defining refinement ratios between levels
                       m_domain,             // THe entire domain without refinement: the base grid
                       alpha=0.0,            // No identity term
                       beta=1.193E9,         // beta=1/(4*pi*G)
                       m_dx)                 // coarsest grid spacing
 */

// Clean up and close up shop
// Output gravpot [phi] and acceleration [-grad(phi)]: used by PLUTO/Src/HD/prim_eqn.c
