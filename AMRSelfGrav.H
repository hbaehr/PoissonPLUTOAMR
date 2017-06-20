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

// The following parts go in PatchPluto.H header file
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
/*
