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

#include "CH_Attach.H"

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

class AMRPoissonPluto : public AMRPoissonOp
{
public:
  /// Constructor
  /**
  */
  AMRPoissonPluto();

  /// Destructor
  /**
  */
  virtual ~AMRPoissonPluto();

  /// Define the parameters the object needs
  /**
  */
  virtual void define(Vector<LevelData<FArrayBox>* >    a_rhs,
                      Vector<DisjointBoxLayout>         a_allGrids,
                      Vector<ProblemDomain>             a_domain,
                      Vector<Real>                      a_dx,
                      Vector<int>                       a_ref_ratio,
                      int                               a_numLevels);

  /*extern void convDiriBC_RBGS( FArrayBox&              a_state,
                               const FArrayBox&        a_rhs,
                               const Box&              a_valid,
                               const ProblemDomain&    a_domain,
                               Real                    a_dx,
                               int                     a_whichpass,
                               int                     a_dir,
                               Side::LoHiSide          a_side);

  /*void convergeGS_BC( FArrayBox& a_state,
                      const FArrayBox& a_rhs,
                      const Box& a_valid,
                      const ProblemDomain& a_domain,
                      Real a_dx,
                      int a_whichpass )

/*  virtual void setRHS(Vector<LevelData<FArrayBox>* > a_rhs,
                      Vector<ProblemDomain>&         a_domain,
                      Vector<int>&                   a_ref_ratio,
                      Vector<Real>&                  a_dx,
                      int                            a_level);*/

  virtual void setupSolver(AMRMultiGrid<LevelData<FArrayBox> >       *a_amrSolver,
                           LinearSolver<LevelData<FArrayBox> >&       a_bottomSolver,
                           Vector<DisjointBoxLayout>&                 a_grids,
                           Vector<ProblemDomain>&                     a_domain,
                           Vector<int>&                               a_ref_ratio,
                           Vector<Real>&                              a_dx,
                           int                                        a_level);

  virtual Vector<LevelData<FArrayBox>* > runSolver();

  Vector<LevelData<FArrayBox>* > phi;

  Vector<LevelData<FArrayBox>* > m_phi;

protected:

  Vector<ProblemDomain> m_domain;

  int m_numLevels;

  Vector<Real> m_dx;

  Vector<int> m_ref_ratio;

  Vector<DisjointBoxLayout> m_grids;

  Vector<LevelData<FArrayBox>* > m_rhs;

  bool m_isDefined;

  int maxLevels;

  int numLevels;
};
