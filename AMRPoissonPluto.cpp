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
//using std::ifstream;
//using std::ios;

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
#include "memtrack.H"

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BCFunc.H"
#include "BiCGStabSolver.H"
#include "CONSTANTS.H"

#include "AMRPoissonOp.H"
#include "AMRMultiGrid.H"

#ifdef CH_LINUX
// Should be undefined by default
#define TRAP_FPE
#undef  TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

// Flag everything as not defined or set
AMRPoissonPluto::AMRPoissonPluto()
{
  m_isDefined = false;
}

AMRPoissonPluto::~AMRPoissonPluto()
{
}
/*
class GlobalBCRS
{
public:
  static std::vector<bool> s_printedThatLo, s_printedThatHi;
  static std::vector<int> s_bcLo, s_bcHi;
  static RealVect s_trigvec;
  static bool s_areBCsParsed, s_valueParsed, s_trigParsed;
};

std::vector<bool> GlobalBCRS::s_printedThatLo = std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi = std::vector<bool>(SpaceDim, false);
std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
RealVect          GlobalBCRS::s_trigvec = RealVect::Zero;
bool              GlobalBCRS::s_areBCsParsed= false;
bool              GlobalBCRS::s_valueParsed= false;
bool              GlobalBCRS::s_trigParsed= false;*/

// Define this object and the boundary condition object
void AMRPoissonPluto::define(Vector<LevelData<FArrayBox>* >   a_rhs,
                             Vector<DisjointBoxLayout>        a_grids,
                             Vector<ProblemDomain>            a_domain,
                             Vector<Real>                     a_dx,
                             Vector<int>                      a_ref_ratio,
                             int                              a_numLevels)
{
  // Store the level data to be used later
  m_rhs       = a_rhs;
  m_grids     = a_grids;
  m_domain    = a_domain;
  m_dx        = a_dx;
  m_numLevels = a_numLevels;
  m_ref_ratio = a_ref_ratio;
  m_isDefined = true;
}

int s_verbosity = 1;

//  -----------------------------------------
// boundary condition stuff
//  -----------------------------------------
///
/**
 */

void ParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  a_values[0]=0.01;
}

void ParseBC(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             bool a_homogeneous)
{

  if (!a_domain.domainBox().contains(a_state.box()))
    {/*
      if (!GlobalBCRS::s_areBCsParsed)
        {
          //ParmParse pp;
          //pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
          //pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
          for (int i=0; i<CH_SPACEDIM; ++i)
            {
              GlobalBCRS::s_bcLo[i] = 0;
              GlobalBCRS::s_bcHi[i] = 0;
            }
          GlobalBCRS::s_areBCsParsed = true;
        }*/
      Box valid = a_valid;
      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))
            {
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if (!a_domain.domainBox().contains(ghostBoxLo))
                {/*
                   if (GlobalBCRS::s_bcLo[i] == 1)
                     {
                       if (!GlobalBCRS::s_printedThatLo[i])
                         {
                           if (a_state.nComp() != 1)
                             {
                               MayDay::Error("using scalar bc function for vector");
                             }
                           GlobalBCRS::s_printedThatLo[i] = true;
                           if (s_verbosity>5)pout() << "const neum bcs lo for direction " << i << endl;
                         }
                       NeumBC(a_state,
                              valid,
                              a_dx,
                              a_homogeneous,
                              ParseValue,
                              i,
                              Side::Lo);
                     }
                   else if (GlobalBCRS::s_bcLo[i] == 0)
                     {
                       if (!GlobalBCRS::s_printedThatLo[i])
                         {
                           if (a_state.nComp() != 1)
                             {
                               MayDay::Error("using scalar bc function for vector");
                             }
                           GlobalBCRS::s_printedThatLo[i] = true;
                           if (s_verbosity>5)pout() << "const diri bcs lo for direction " << i << endl;
                         }*/
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             true,
                             ParseValue,
                             i,
                             Side::Lo,
                             1);
                    /*}
                  else
                    {
                      MayDay::Error("bogus bc flag lo");
                    }*/
                }

              if (!a_domain.domainBox().contains(ghostBoxHi))
                {/*
                   if (GlobalBCRS::s_bcHi[i] == 1)
                     {
                       if (!GlobalBCRS::s_printedThatHi[i])
                         {
                           if (a_state.nComp() != 1)
                             {
                               MayDay::Error("using scalar bc function for vector");
                             }
                           GlobalBCRS::s_printedThatHi[i] = true;
                           if (s_verbosity>5)pout() << "const neum bcs hi for direction " << i << endl;
                         }
                       NeumBC(a_state,
                              valid,
                              a_dx,
                              a_homogeneous,
                              ParseValue,
                              i,
                              Side::Hi);
                     }
                   else if (GlobalBCRS::s_bcHi[i] == 0)
                     {
                       if (!GlobalBCRS::s_printedThatHi[i])
                         {
                           if (a_state.nComp() != 1)
                             {
                               MayDay::Error("using scalar bc function for vector");
                             }
                           GlobalBCRS::s_printedThatHi[i] = true;
                           if (s_verbosity>5)pout() << "const diri bcs hi for direction " << i << endl;
                         }*/
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             true,
                             ParseValue,
                             i,
                             Side::Hi,
                             1);
                  /*  }
                  else
                    {
                       MayDay::Error("bogus bc flag hi");
                    }*/
                }
            } // end if is not periodic in ith direction
        }
    }
}
#if 0
extern void convDiriBC_RBGS( FArrayBox&      a_state,
                             const FArrayBox& a_rhs,
                             const Box&      a_valid,
                             const ProblemDomain& a_domain,
                             Real a_dx,
                             int             a_whichpass,
                             int             a_dir,
                             Side::LoHiSide  a_side);

void convergeGS_BC( FArrayBox& a_state,
                    const FArrayBox& a_rhs,
                    const Box& a_valid,
                    const ProblemDomain& a_domain,
                    Real a_dx,
                    int a_whichpass )
{

  if (!a_domain.domainBox().contains(a_state.box()))
    {
      if (!GlobalBCRS::s_areBCsParsed)
        {
          //ParmParse pp;
          //pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
          //pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
          for (int i=0; i<CH_SPACEDIM; ++i)
            {
              GlobalBCRS::s_bcLo[i] = 0;
              GlobalBCRS::s_bcHi[i] = 0;
            }
          GlobalBCRS::s_areBCsParsed = true;
        }

      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))
            {
              Box ghostBoxLo = adjCellBox(a_valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(a_valid, i, Side::Hi, 1);
              if ( !a_domain.domainBox().contains(ghostBoxLo) && GlobalBCRS::s_bcLo[i]==0 )
                {
                  Interval stateInterval = a_state.interval();
                  // iterate
                  for (int kk=0;kk<10;kk++)
                    {
                      // apply BC
                      DiriBC(a_state,
                             a_valid,
                             a_dx,
                             true,
                             ParseValue,
                             i,
                             Side::Lo,
                             1);
                      // apply G-S to boundary
                      convDiriBC_RBGS(a_state,
                                      a_rhs,
                                      a_valid,
                                      a_domain,
                                      a_dx,
                                      a_whichpass,
                                      i,
                                      Side::Lo);
                    }
                }
              if (!a_domain.domainBox().contains(ghostBoxHi) && GlobalBCRS::s_bcHi[i]==0 )
                {
                  Interval stateInterval = a_state.interval();
                  // iterate
                  for (int kk=0;kk<10;kk++)
                    {
                      // apply BC
                      DiriBC(a_state,
                             a_valid,
                             a_dx,
                             true,
                             ParseValue,
                             i,
                             Side::Hi,
                             1);
                      // apply G-S to boundary
                      convDiriBC_RBGS(a_state,
                                      a_rhs,
                                      a_valid,
                                      a_domain,
                                      a_dx,
                                      a_whichpass,
                                      i,
                                      Side::Hi);
                    }
                }
            } // end if is not periodic in ith direction
        }
    }
}
#endif
void AMRPoissonPluto::setupSolver(AMRMultiGrid<LevelData<FArrayBox> > *a_amrPoissonSolver,
                                  LinearSolver<LevelData<FArrayBox> >& a_bottomSolver,
                                  Vector<DisjointBoxLayout>&           a_grids,
                                  Vector<ProblemDomain>&               a_domain,
                                  Vector<int>&                         a_ref_ratio,
                                  Vector<Real>&                        a_dx,
                                  int                                  a_level)
{
   CH_TIME("setupSolver");

   AMRPoissonOpFactory opFactory;

   // solving poisson problem here
   Real alpha = 0.0;
   Real beta = 1.0/(4*3.14159265);

   opFactory.define(a_domain[0],
                    a_grids,
                    a_ref_ratio,
                    a_dx[0],
                    &ParseBC, alpha, beta);

   AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;

   a_amrPoissonSolver->define(a_domain[0], castFact,
                      &a_bottomSolver, a_level);

   // Multigrid solver parameters
   int numSmooth = 4;
   int numMG = 1;                                                                // Type of multigrid solver: 1=V-cycle, 2=W-cycle, etc.
   int maxIter = 100;
   Real eps = 1.0e-9;
   Real hang = 1.0e-10;

   Real normThresh = 1.0e-30;
   a_amrPoissonSolver->setSolverParameters(numSmooth, numSmooth, numSmooth,
                                numMG, maxIter, eps, hang, normThresh);
 }

int AMRPoissonPluto::runSolver()
{
   CH_TIME("runSolver");

   int status = 0, mg_type = 0;

   //int s_verbosity = 4;

   // allocate solution and RHS, initialize RHS
   int numLevels = grids.size();
   Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);
   //Vector<LevelData<FArrayBox>* > rhs(numLevels, NULL);
   // this is for convenience
   Vector<LevelData<FArrayBox>* > resid(numLevels, NULL);

   for (int lev=0; lev<=numLevels-1; lev++)
     {
       const DisjointBoxLayout& levelGrids = m_grids[lev];
       phi[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);
       //rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
       resid[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
     }

   // initialize solver
   AMRMultiGrid<LevelData<FArrayBox> > *amrPoissonSolver;
   if ( mg_type==0 )
     {
       amrPoissonSolver = new AMRMultiGrid<LevelData<FArrayBox> >();
     }
   else
     {
       MayDay::Error("FAS not supported");
     }

   BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
   //bottomSolver.m_verbosity = s_verbosity-2;
   setupSolver(amrPoissonSolver, bottomSolver, m_grids, m_domain,
               m_ref_ratio, m_dx, m_numLevels);

   // do solve
   int iterations = 1;

   for (int iiter = 0; iiter < iterations; iiter++)
     {
       bool zeroInitialGuess = true;
       pout() << "about to go into solve" << endl;
       amrPoissonSolver->solve(phi, m_rhs, numLevels-1, 0, zeroInitialGuess);
       pout() << "done solve" << endl;
     }

   // write results to file
   bool writePlots = false;
//   ppMain.query("writePlotFiles", writePlots);

#ifdef CH_USE_HDF5

   if (writePlots)
     {
       int numLevels = m_grids.size();
       Vector<LevelData<FArrayBox>* > plotData(numLevels, NULL);

       pout() << "Write Plots. norm=" << amrPoissonSolver->computeAMRResidual(resid,phi,m_rhs,numLevels-1,0) << endl;

       for (int lev=0; lev<numLevels; lev++)
         {
           plotData[lev] = new LevelData<FArrayBox>(amrGrids[lev],
                                                    3, IntVect::Zero);

           Interval phiInterval(0,0);
           phi[lev]->copyTo(phiInterval, *plotData[lev], phiInterval);
           Interval rhsInterval(1,1);
           m_rhs[lev]->copyTo(phiInterval, *plotData[lev], rhsInterval);
           Interval resInterval(2,2);
           resid[lev]->copyTo(phiInterval, *plotData[lev], resInterval);
         }

       string fname = "poissonOut.";

       char suffix[30];
       sprintf(suffix, "%dd.hdf5",SpaceDim);
       fname += suffix;

       Vector<string> varNames(3);
       varNames[0] = "phi";
       varNames[1] = "rhs";
       varNames[2] = "res";

       Real bogusVal = 1.0;

       WriteAMRHierarchyHDF5(fname,
                             m_grids,
                             plotData,
                             varNames,
                             m_domain[0].domainBox(),
                             m_dx[0],
                             bogusVal,
                             bogusVal,
                             m_ref_ratio,
                             numLevels);

       // clean up
       for (int lev=0; lev<plotData.size(); lev++)
         {
           delete plotData[lev];
         }
     } // end if writing plots
#endif // end if HDF5

   delete amrPoissonSolver;

   return status;

   // clean up
   for (int lev=0; lev<phi.size(); lev++)
     {
       delete phi[lev];
       //delete m_rhs[lev];
       delete resid[lev];
     }
}
