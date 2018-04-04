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

// Other things which need to be addressed:
// 1) Boundary conditions
// 2) Fill the ghost zones of m_gravpot or a_gravpot by interpolation? see LevelPluto for more info
// 3) Go through each input from PLUTO and decide whether they should be  temporary 'm_' or permanent 'a_' for now they are all 'a_'

/*
* A few notes:
* m_ indicates a temporary storage structure and should be used within routines
*    and re-used each timestep
*
* a_ indicates the main which stores data between timesteps
*
*
*
*/// Flag everything as not defined or set
AMRPoissonPluto::AMRPoissonPluto()
{
  m_isDefined = false;
  m_isBoundarySet = false;
  m_isRiemannSet = false;
  m_isCurrentTimeSet = false;
  m_isCurrentBoxSet = false;
}

AMRPoissonPluto::~AMRPoissonPluto()
{
}

// Define this object and the boundary condition object
void AMRPoissonPluto::define(ProblemDomain& a_domain,
                          const Real&    a_dx,
                          const int&     a_level,
                          const int&     a_numGhost)
{

  // Store the domain and grid spacing
  m_domain = a_domain;
  m_dx = a_dx;
  m_level = a_level;
  m_numGhost = a_numGhost;
  m_isDefined = true;
}

int s_verbosity = 1;

//  -----------------------------------------
// boundary condition stuff                                                      // Are these boundary conditions between AMR levels or for the entire domain?
//  -----------------------------------------
///
/**
 */
class GlobalBCRS                                                                 // Is there any documentation on GlobalBCRS?
{
public:
  static std::vector<bool> s_printedThatLo, s_printedThatHi;                     // why bool?
  static std::vector<int> s_bcLo, s_bcHi;                                        // Lo and Hi refer to what exactly?
  static RealVect s_trigvec;
  static bool s_areBCsParsed, s_valueParsed, s_trigParsed;
};

std::vector<bool> GlobalBCRS::s_printedThatLo = std::vector<bool>(SpaceDim, false); // vector<bool> is ???
std::vector<bool> GlobalBCRS::s_printedThatHi = std::vector<bool>(SpaceDim, false); //
std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
RealVect          GlobalBCRS::s_trigvec = RealVect::Zero;
bool              GlobalBCRS::s_areBCsParsed= false;
bool              GlobalBCRS::s_valueParsed= false;
bool              GlobalBCRS::s_trigParsed= false;
Real bcVal=0.0;

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

 void AMRLevelPluto::setupSolver(AMRMultiGrid<LevelData<FArrayBox> > *a_amrPoissonSolver,              // Name of the solver
                                 LinearSolver<LevelData<FArrayBox> >& a_bottomSolver,                // The bottom solver is what solves the Poisson equation on the coarsest level
                                 Vector<DisjointBoxLayout>& a_grids,                           // Grids for each AMR level (i.e. there are multiple grids per level)
                                 Vector<ProblemDomain>& a_domain,                              // Entire domain
                                 Vector<int>& a_ref_ratio,                                     // Refinement ratios between levels
                                 Vector<Real>& a_dx,                                           // grid spacing
                                 int a_level)                                                        // number of most refined level
 {
   AMRPoissonOpFactory opFactory;                                                // Create an instance of AMRPoissonOpFactory

   // solving poisson problem here
   Real alpha = 0.0;                                                             // Constants which determine form of the Poisson equation
   Real beta = 1.0/(4*3.14159265);                                               // \beta will be 1/(4*\pi*G) in cgs or code units? For now, G=1

   opFactory.define(a_domain[0],                                                 // Define the parameters that go into each instance of opFactory
                    a_grids,                                                     // Defining these parameters is the whole point of setupGrids, but these
                    a_ref_ratio,                                                 // are already defined by the main part of PLUTO-Chombo
                    a_dx[0],                                                     // so setupGrids is probably unnecessary.
                    &ParseBC, alpha, beta);

   AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;

   a_amrPoissonSolver->define(a_domain[0], castFact,                                    // Define
                      &a_bottomSolver, numLevels);

   // Multigrid solver parameters
   int numSmooth = 4;
   int numMG = 1;                                                                // Type of multigrid solver: 1=V-cycle, 2=W-cycle, etc.
   int maxIter = 100;
   Real eps = 1.0e-9;
   Real hang = 1.0e-10;

   Real normThresh = 1.0e-30;                                                    // What is this threshold? Look to the a_amrSolver
   a_amrPoissonSolver->setSolverParameters(numSmooth, numSmooth, numSmooth,             // Input all the parameters
                                numMG, maxIter, eps, hang, normThresh);
 }

 int runSolver()                                                                 // Now that everything is set up calculate the potential
 {
   int status = 0, mg_type = 0;                                                  // ???

   //int s_verbosity = 4;

   Vector<DisjointBoxLayout> grids;
   Vector<ProblemDomain> domain;
   Vector<int> ref_ratio;
   Vector<Real> dx;
   int level;

   // allocate solution and RHS, initialize RHS
   int numLevels = grids.size();                                                 // Set up containers
   Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);                          // \phi container
   Vector<LevelData<FArrayBox>* > rhs(numLevels, NULL);                          // \rho container
   // this is for convenience
   Vector<LevelData<FArrayBox>* > resid(numLevels, NULL);                        // residual container

   for (int lev=0; lev<=maxLevel; lev++)                                            // Not sure how to loop over the AMR levels when PLUTO doesn't, maybe a different container?
     {
       const DisjointBoxLayout& levelGrids = grids[lev];                         // lev = level dummy
       phi[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);        // create space for \phi data for each level
       rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);        // create space for \rho date for each level
       resid[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);      // create space for residual for each level
     }

   setRHS(rhs, domain, ref_ratio, dx, level);                                   // set the RHS, wasn't this done in the beginning already?

   // initialize solver
   AMRMultiGrid<LevelData<FArrayBox> > *amrPoissonSolver;                               // Not a_amrSolver? Where is amrSolver defined? And why a pointer?
   if ( mg_type==0 )
     {
       amrSolver = new AMRMultiGrid<LevelData<FArrayBox> >();                    // This is the solver from one level to another
     }
   else
     {
       MayDay::Error("FAS not supported");
     }

   BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;                           // So a v-cycle needs a direct solver on the coarsest level
   //bottomSolver.m_verbosity = s_verbosity-2;
   setupSolver(amrPoissonSolver, bottomSolver, amrGrids, amrDomains,
               refRatios, amrDx, level);
                                                                                 // or was that just the rules and parameters and this is the execution?
   // do solve
   int iterations = 1;

   for (int iiter = 0; iiter < iterations; iiter++)                              // interate over
     {
       bool zeroInitialGuess = true;
       pout() << "about to go into solve" << endl;
       amrPoissonSolver->solve(phi, rhs, level, 0, zeroInitialGuess);           // Here is where it is all put together
       pout() << "done solve" << endl;
     }

   // write results to file

/*
*  This may be useful for testing an debugging purposes, but in the end,
*  I would prefer to keep all the information in one place with the rest of
*  the conserved/primative variables of PLUTO
*
*/

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
