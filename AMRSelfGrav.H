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

class AMRPoissonPluto : public AMRPoisson
{
public:
  /// Constructor
  /**
  */
  AMRPoissonPluto();

  /// Destructor
  /**
  */
  virtual ~AMRLevelPluto();

  /// Define the parameters the object needs
  /**
  */
  void defineParams(
  );

virtual void ParseValue(Real* pos,                                                       // Looks like this would normally parse a BC value settings from a file
                        int* dir,                                                        // dir = direction? as in dimension 1,2,3?
                        Side::LoHiSide* side,                                            // Which side of the coordinate boundary?
                        Real* a_values);

virtual void ParseBC(FArrayBox& a_state,                                                 // See documentation for BCfunction. State is the value of Phi within box
                     const Box& a_valid,                                                 // domain of the box (might be a subdomain?) around which the BCs are 'built'
                     const ProblemDomain& a_domain,                                      // This is the coarsest domain
                     Real a_dx,                                                          // without grid data is this just the coarse dx
                     bool a_homogeneous);

virtual void setRHS(Vector<LevelData<FArrayBox>* > rhs,                                  // Output array of \rho: a_U[RHO]?
                    Vector<ProblemDomain>& a_domain,                                     // Grid domain
                    Vector<int>& a_ref_ratio,                                            // Refinement ratios between levels
                    Vector<Real>& a_dx,                                                  // dx: grid spacing
                    int a_level);                                                        // number of most refined level = m_level

virtual int runSolver();

virtual void setupSolver(AMRMultiGrid<LevelData<FArrayBox> >        *a_amrSolver,           // Name of the Multigrid solver
                         LinearSolver<LevelData<FArrayBox> >&       a_bottomSolver,         // Solver at the bottom of the V-cycle
                         Vector<DisjointBoxLayout>&                 a_grids,                // Grids for each AMR level
                         Vector<ProblemDomain>&                     a_domain,               // Entire domain
                         Vector<int>&                               a_ref_ratio,            // Refinement ratios between levels
                         Vector<Real>&                              a_dx,                // grid spacing
                         int                                        a_level);         // number of most refined level
