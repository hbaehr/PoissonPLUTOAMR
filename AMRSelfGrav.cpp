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

#include "UsingNamespace.H"

#ifdef CH_LINUX
// Should be undefined by default
#define TRAP_FPE
#undef  TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

/* List of things which need to happen
 *
 * Calculate self-gravitational potential
 * Calculate aceleration due to self-gravity
 * Export calculations to main part of PLUTO
 *  - Store calculations somewhere (temporarily)?
 * 
 *
 */
    
// Use CHOMBO functions to solve for the self-gravitational potential
    
// Define a host of PoissonOps for all necessary levels
void AMRPoissonOpFactory::define(a_coarseDomain   // domain at the coarsest level
                                 a_grids          // AMR heirarchy
                                 a_refRatios      // refinement ratios between levels
                                 a_coarsedx       // grid spacing at the coarsest level
                                 a_bc             // boundary conditions
                                 a_alpha=0.0      // identity coefficient = 0
                                 a_beta=1.0       // Lapacian coefficient = 1/(4*pi*G)
                            );

//
void AMRPoissonOp::(
                    );

//
virtual void AMRPoissonOp::applyOp(
                                   );

// Getting an accurate value
virtual void AMRPoissonOp::residual(LevelData<FArrayBox>&  	a_lhs,      // left-hand side
                                    const LevelData<FArrayBox>&  a_phi, // self-gravitational potential
                                    const LevelData<FArrayBox>&  a_rhs, // right-hand side (rho)
                                    bool  a_homogeneous = false         // 
                                    );

// Example from the documentation

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

void solveElliptic(Vector<LevelData<FArrayBox>* >& phi,
                   const Vector<LevelData<FArrayBox>* > rhs,
                   const Vector<DisjointBoxLayout>& grids,
                   const Vector<int>& refRatios,
                   const ProblemDomain& level0Domain,
                   Real alpha, Real beta, Real coarsestDx)
{
int numlevels = rhs.size();

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

// Calculate physics of self-gravity grad(phi) = g (Do this here or in the main PLUTO code?)
// It might be useful to do it here because Chombo should have no issues doing the gradient 
// over grid refinement levels

// Export calculations to the main part of PLUTO (phi and grad(phi)?)
 
// Fluxes?

// Clean up and close up shop