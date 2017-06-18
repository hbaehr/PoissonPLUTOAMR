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
    CH_TIME("setup::Udefine");
    m_gravpot.define(m_grids,m_numCons,m_numGhost*IntVect::Unit);
  }

// And maybe fill the ghost zones of m_gravpot or a_gravpot by interpolation? see LevelPluto for more info

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
