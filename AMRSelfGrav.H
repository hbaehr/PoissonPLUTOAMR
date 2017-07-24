/*
 *   Header file for speeding up compilation
 *   This header is only necessary if AMRSelfgrav is a class rather than a function
 */


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

virtual void initiate(LevelData<FArrayBox>& a_gravpot);

virtual void starter(LevelData<FArrayBox>& a_gravpot);

virtual void solveSelfGravPot(Vector<LevelData<FArrayBox>*>&       a_gravpot,
                              Vector<LevelData<FArrayBox>*>        a_U,
                              const Vector<DisjointBoxLayout>&     a_grids,
                              const Vector<int>&                   m_ref_ratio,
                              const ProblemDomain&                 a_domain,
                              Real                                 alpha=0.0,
                              Real                                 beta=1.193E9, // This should perhaps be problem defined
                              Real                                 a_dx);

/* **********************************************************
void solveElliptic(Vector<LevelData<FArrayBox>* >& phi,         // Output self-gravity potential: m_gravpot
                   const Vector<LevelData<FArrayBox>* > rhs,    // Input density: m_UNew
                   const Vector<DisjointBoxLayout>& grids,      // Grid geometries: m_ref_ratio
                   const Vector<int>& refRatios,                // Vector defining refinement ratios between levels
                   const ProblemDomain& level0Domain,           //
                   Real alpha, Real beta, Real coarsestDx)      // constants alpha=0, beta=1/(4*pi*G)
*/

virtual int runSolver();

virtual void setupSolver(AMRMultiGrid<LevelData<FArrayBox> >        *a_amrSolver,           // Name of the solver
                         LinearSolver<LevelData<FArrayBox> >&       a_bottomSolver,         //
                         const Vector<DisjointBoxLayout>&           a_grids,                // Grids for each AMR level
                         const Vector<ProblemDomain>&               a_domain,               // Entire domain
                         const Vector<int>&                         a_ref_ratio,            // Refinement ratios between levels
                         const Vector<Real>&                        a_amrDx,                // *** dx: not sure what this value is atm
                         int                                        a_finestLevel);         // *** number of most refined level

virtual void setRHS(Vector<LevelData<FArrayBox>* >                  a_rhs,                  // Output array of \rho
                    Vector<ProblemDomain>&                          a_domain,               // Grid domain
                    Vector<int>&                                    m_ref_ratio,            // Refinement ratios between levels
                    Vector<Real>&                                   a_amrDx,                // *** dx: not sure what this value is atm
                    int                                             a_finestLevel)          // *** number of most refined level
