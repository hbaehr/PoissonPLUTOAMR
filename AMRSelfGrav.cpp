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
    
// Example from the Chombo documentation

/*
 * solveElliptic solves (alpha I + beta Laplacian) phi = rhs
 * using AMRMultiGrid and AMRPoissonOp
 * Inputs:
 * rhs: Right-hand side of the solve over the level.
 * grids: AMRHierarchy of grids
 * refRatio: refinement ratios
 * level0Domain: domain at the coarsest AMR level
 * coarsestDx: grid spacing at the coarsest level
 * alpha: identity coefficient
 * beta: Laplacian coefficient
 * Outputs:
 * phi = (alpha I + beta Lapl)^{-1}(rhs)
 */

void solveElliptic(Vector<LevelData<FArrayBox>* >& phi,         // Output self-gravity potential
                   const Vector<LevelData<FArrayBox>* > rhs,    // Input density 
                   const Vector<DisjointBoxLayout>& grids,      // Grid geometries
                   const Vector<int>& refRatios,                // Vector defining refinement ratios between levels
                   const ProblemDomain& level0Domain,           // 
                   Real alpha, Real beta, Real coarsestDx)      // constants alpha=0, beta=1/(4*pi*G)
{
int numlevels = rhs.size(); // A different array for each refinement level

//define the operator factory
AMRPoissonOpFactory opFactory;
opFactory.define(level0Domain,
                 grids, refRatios, coarsestDx,
                 &ParseBC, alpha, beta);

//this is the solver we shall use
AMRMultiGrid<LevelData<FArrayBox> > solver;

//this is the solver for the bottom of the muligrid v-cycle
BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;

//bicgstab can be noisy
bottomSolver.m_verbosity = 0;

//define the solver
solver.define(level0Domain, opFactory, &bottomSolver, numlevels);

//we want to solve over the whole hierarchy
int lbase = 0;

//so solve already.
solver.solve(phi, rhs, numlevels-1, lbase);
}
 
// Fluxes?

// Clean up and close up shop