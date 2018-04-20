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

// Define this object and the boundary condition object
void AMRPoissonPluto::defineParams(ProblemDomain&      a_domain,
                                   Real&               a_dx,
                                   int                 a_level,
                                   int&                a_ref_ratio)
{

  // Store the domain and grid spacing
  m_domain = a_domain;
  m_dx = a_dx;
  m_level = a_level;
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

void AMRPoissonPluto::ParseValue(Real* pos,
                                 int* dir,
                                 Side::LoHiSide* side,
                                 Real* a_values)
{
  a_values[0]=0.;
}

void AMRPoissonPLuto::ParseBC(FArrayBox& a_state,
                              const Box& a_valid,
                              const ProblemDomain& a_domain,
                              Real a_dx,
                              bool a_homogeneous)
{

  if (!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))
            {
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if (!a_domain.domainBox().contains(ghostBoxLo))
                {
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             true,
                             ParseValue,
                             i,
                             Side::Lo,
                             1);
                }

              if (!a_domain.domainBox().contains(ghostBoxHi))
                {
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             true,
                             ParseValue,
                             i,
                             Side::Hi,
                             1);
                }
            } // end if is not periodic in ith direction
        }
    }
}

void AMRPoissonPluto::setRHS(Vector<LevelData<FArrayBox>* > a_rhs,
                             Vector<ProblemDomain>&         a_domain,
                             Vector<int>&                   a_ref_ratio,
                             Vector<Real>&                  a_dx,
                             int                            a_level)
{
   int maxLevel = 0;
   maxLevel = atoi(ParamFileGet("Levels",1));
   int numLevels = MAX(maxLevel, 1)

   // Innitialize rhs container for Poisson solver
   a_rhs.resize(maxLevel);

   // Setup variable for the conserved variable rho to copy it into the new structure
   Interval rhsInterval(0,0);
   m_UOld.copyTo(rhsInterval,*a_rhs[a_level],rhsInterval);

   // Resize Vector<> containers for Poisson rhs
   a_grids.resize(maxLevel);
   a_domain.resize(maxLevel);
   a_ref_ratio.resize(maxLevel);
   a_dx.resize(maxLevel);

   // Put values for current level in proper container for the solver
   a_grids[a_level] = m_grids;
   a_domain[a_level] = m_problem_domain;
   a_ref_ratio[a_level] = m_ref_ratio;
   a_dx[a_level] = m_dx;
 }

 void AMRLevelPluto::setupSolver(AMRMultiGrid<LevelData<FArrayBox> > *a_amrPoissonSolver,
                                 LinearSolver<LevelData<FArrayBox> >& a_bottomSolver,
                                 Vector<DisjointBoxLayout>&           a_grids,
                                 Vector<ProblemDomain>&               a_domain,
                                 Vector<int>&                         a_ref_ratio,
                                 Vector<Real>&                        a_dx,
                                 int                                  a_level)
 {
   AMRPoissonOpFactory opFactory;                                                // Create an instance of AMRPoissonOpFactory

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
                      &a_bottomSolver, numLevels);

   // Multigrid solver parameters
   int numSmooth = 4;
   int numMG = 1;                                                                // Type of multigrid solver: 1=V-cycle, 2=W-cycle, etc.
   int maxIter = 100;
   Real eps = 1.0e-9;
   Real hang = 1.0e-10;

   Real normThresh = 1.0e-30;                                                    // What is this threshold? Look to AMRMultiGrid
   a_amrPoissonSolver->setSolverParameters(numSmooth, numSmooth, numSmooth,
                                numMG, maxIter, eps, hang, normThresh);
 }

 int runSolver()
 {
   int status = 0, mg_type = 0;

   //int s_verbosity = 4;

   Vector<DisjointBoxLayout> grids;
   Vector<ProblemDomain> domain;
   Vector<int> ref_ratio;
   Vector<Real> dx;
   int level;

   // allocate solution and RHS, initialize RHS
   int numLevels = grids.size();
   Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);
   Vector<LevelData<FArrayBox>* > rhs(numLevels, NULL);
   // this is for convenience
   Vector<LevelData<FArrayBox>* > resid(numLevels, NULL);

   for (int lev=0; lev<=maxLevel; lev++)
     {
       const DisjointBoxLayout& levelGrids = grids[lev];
       phi[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);
       rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
       resid[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
     }

   // initialize solver
   AMRMultiGrid<LevelData<FArrayBox> > *amrPoissonSolver;
   if ( mg_type==0 )
     {
       amrSolver = new AMRMultiGrid<LevelData<FArrayBox> >();
     }
   else
     {
       MayDay::Error("FAS not supported");
     }

   BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
   //bottomSolver.m_verbosity = s_verbosity-2;
   setupSolver(amrPoissonSolver, bottomSolver, grids, domain,
               ref_ratio, dx, level);

   setRHS(rhs, domain, ref_ratio, dx, level);

   // do solve
   int iterations = 1;

   for (int iiter = 0; iiter < iterations; iiter++)
     {
       bool zeroInitialGuess = true;
       pout() << "about to go into solve" << endl;
       amrPoissonSolver->solve(phi, rhs, level, 0, zeroInitialGuess);
       pout() << "done solve" << endl;
     }

   // write results to file
   bool writePlots = false;
//   ppMain.query("writePlotFiles", writePlots);

#ifdef CH_USE_HDF5

   if (writePlots)
     {
       int numLevels = level +1;
       Vector<LevelData<FArrayBox>* > plotData(numLevels, NULL);

       pout() << "Write Plots. norm=" << amrPoissonSolver->computeAMRResidual(resid,phi,rhs,level,0) << endl;

       for (int lev=0; lev<numLevels; lev++)
         {
           plotData[lev] = new LevelData<FArrayBox>(amrGrids[lev],
                                                    3, IntVect::Zero);

           Interval phiInterval(0,0);
           phi[lev]->copyTo(phiInterval, *plotData[lev], phiInterval);
           Interval rhsInterval(1,1);
           rhs[lev]->copyTo(phiInterval, *plotData[lev], rhsInterval);
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
                             amrGrids,
                             plotData,
                             varNames,
                             amrDomains[0].domainBox(),
                             amrDx[0],
                             bogusVal,
                             bogusVal,
                             refRatios,
                             numLevels);

       // clean up
       for (int lev=0; lev<plotData.size(); lev++)
         {
           delete plotData[lev];
         }
     } // end if writing plots
#endif // end if HDF5

   // clean up
   for (int lev=0; lev<phi.size(); lev++)
     {
       delete phi[lev];
       delete rhs[lev];
       delete resid[lev];
     }

   delete amrPoissonSolver;

   return status;
}

// Output gravpot [phi] and acceleration [-grad(phi)]: used by PLUTO/Src/HD/prim_eqn.c
