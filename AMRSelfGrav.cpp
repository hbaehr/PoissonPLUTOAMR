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

/* For time being, unnecessary
  // initialize some temporary storage for the gravitational potential
  {
  //  CH_TIME("setup::Udefine"); // This is used to time processes for debugging purposes
    m_gravpot.define(m_grids,m_numCons,m_numGhost*IntVect::Unit);
  }
  */

// Other things which need to be addressed:
// 1) Boundary conditions
// 2) Initial conditions
// 2) Fill the ghost zones of m_gravpot or a_gravpot by interpolation? see LevelPluto for more info

/* Set the right hand side of Poisson equation
*  \nabla \Phi = 4*\pi*G*\rho
*
*  Since everything on the RHS outside of \rho is constant, these are put into
*  the constant value b. However, since \rho has already been calculate in the
*  previous step, it might be unnecessary to do this step.
*
*  I will include it for now for completeness and perhaps testing purposes
*/
void setRHS(Vector<LevelData<FArrayBox>* > a_rhs,                                // Output array of \rho
            Vector<ProblemDomain>& a_domain,                                     // Grid domain
            Vector<int>& m_ref_ratio,                                            // Refinement ratios between levels
            Vector<Real>& a_dx,                                                  // *** dx: not sure what this value is atm
            int a_finestLevel)                                                   // *** number of most refined level
{
  CH_TIME("setRHS");

  for (int lev=0; lev<=a_finestLevel; lev++)                                     // Looping over all levels
    {
      LevelData<FArrayBox>& levelRhs = *(a_rhs[lev]);                            // Each a_rhs[lev] is an FArrayBox temporarily stored as levelRhs
      const DisjointBoxLayout& levelGrids = levelRhs.getBoxes();                 // Get box indeces for each level

      // rhs is cell-centered...
      RealVect ccOffset = 0.5*a_amrDx[lev]*RealVect::Unit;                       // Adjust for the cell-centered location of the data ???

      DataIterator levelDit = levelGrids.dataIterator();                         // Write data to a level
      for (levelDit.begin(); levelDit.ok(); ++levelDit)                          // Iterate over all points on the level?
        {
          FArrayBox& thisRhs = levelRhs[levelDit];                               // Dummy thisRhs for the iteration of this loop

          if (s_probtype == zeroRHS)                                             // Various patterns for the distribution of RHS
            {
              thisRhs.setVal(0.0);                                               // 0 everywhere
            }
          else if (s_probtype == unityRHS)
            {
              thisRhs.setVal(1.0);                                               // 1 everywhere
            }
          else if (s_probtype == sinusoidal)
            {

              BoxIterator bit(thisRhs.box());                                    // Sine wave
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  RealVect loc(iv);
                  loc *= a_amrDx[lev];
                  loc += ccOffset;
                  loc *= 2.0*Pi;

                  thisRhs(iv,0) = D_TERM(sin(loc[0]),
                                         *cos(loc[1]),
                                         *sin(loc[2]));

                 }

            }
          else if (s_probtype == gaussians)
            {
              int numGaussians = 3;                                              // Gaussian distributions
              Vector<RealVect> center(numGaussians,RealVect::Zero);
              Vector<Real> scale(numGaussians, 1.0);
              Vector<Real> strength(numGaussians, 1.0);

              for (int n=0; n<numGaussians; n++)
                {
                  if (n==0)
                    {
                      strength[0] = 1.0;
                      scale[0] = 1.0e-2;
                      center[0] = 0.25*RealVect::Unit;
                    }
                  else if (n == 1)
                    {
                      strength[1] = 3.0;
                      scale[1] = 1.0e-2;
                      center[1] = RealVect(D_DECL(0.5,0.75, 0.75));
                    }
                  else if (n == 2)
                    {
                      strength[2] = 2.0;
                      scale[2] = 1.0e-2;
                      center[2] = RealVect(D_DECL(0.75,0.5, 0.5));
                    }
                  else
                    {
                      MayDay::Error("too many Gaussian sources attempted");
                    }
                }

              thisRhs.setVal(0.0);

              BoxIterator bit(thisRhs.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  RealVect loc(iv);
                  loc *= a_amrDx[lev];
                  loc += ccOffset;

                  for (int n=0; n<numGaussians; n++)
                    {
                      RealVect dist = loc - center[n];
                      Real radSqr = D_TERM(dist[0]*dist[0],
                                           +dist[1]*dist[1],
                                            +dist[2]*dist[2]);

                       Real val = strength[n]*exp(-radSqr/scale[n]);
                       thisRhs(iv,0) += val;
                     }
                 }
             }
           else
             {
               MayDay::Error("undefined problem type");
             }
         } // end loop over grids on this level
     } // end loop over levels
 }


/*  Setting up the parameters of the solver
*
*
*
*
*/


 void setupSolver(AMRMultiGrid<LevelData<FArrayBox> > *a_amrSolver,              // Name of the solver
             LinearSolver<LevelData<FArrayBox> >& a_bottomSolver,                //
             const Vector<DisjointBoxLayout>& a_grids,                           // Grids for each AMR level
             const Vector<ProblemDomain>& a_domain,                              // Entire domain
             const Vector<int>& a_ref_ratio,                                     // Refinement ratios between levels
             const Vector<Real>& a_dx,                                           // *** dx: not sure what this value is atm
             int a_finestLevel)                                                  // *** number of most refined level
 {
   CH_TIME("setupSolver");                                                       // Timing diagnostic

   ParmParse ppSolver("solver");                                                 // ??? WHere does ppSolver come from? Parse parameter information from input files

   int numLevels = a_finestLevel+1;                                              // 0 index count of total levels

   AMRPoissonOpFactory opFactory;                                                // Create an instance of AMRPoissonOpFactory

   // solving poisson problem here
   Real alpha =0.0;                                                              // Constants which determine form of the Poisson equation
   Real beta = 1.0;                                                              // \beta will be 1/(4*\pi*G) in cgs or code units?

   opFactory.define(a_domain[0],                                                 // Define the parameters that go into each instance of opFactory
                    a_grids,
                    a_ref ratio,
                    a_amrDx[0],
                    &ParseBC, alpha, beta);

   AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;  // ??? dummy?

   a_amrSolver->define(a_domain[0], castFact,                                    // Define
                      &a_bottomSolver, numLevels);

   // multigrid solver parameters                                                // Parameters for solving over multiple levels ???
   int numSmooth, numMG, maxIter;                                                // So ParmParse look into the 'inputs' file and for in the instance
   Real eps, hang;                                                               // ppSolver, extracts all values on lines which start with 'solver.'
   ppSolver.get("num_smooth", numSmooth);                                        // and returns the values for the respective variable names
   ppSolver.get("num_mg",     numMG);
   ppSolver.get("max_iterations", maxIter);
   ppSolver.get("tolerance", eps);
   ppSolver.get("hang",      hang);

   Real normThresh = 1.0e-30;
   a_amrSolver->setSolverParameters(numSmooth, numSmooth, numSmooth,
                                numMG, maxIter, eps, hang, normThresh);
   a_amrSolver->m_verbosity = s_verbosity-1;

   // optional parameters
   ppSolver.query("num_pre", a_amrSolver->m_pre);
   ppSolver.query("num_post", a_amrSolver->m_post);
   ppSolver.query("num_bottom", a_amrSolver->m_bottom);
 }

 int runSolver()                                                                 // Now that everything is set up (???) calculate the potential
 {
   CH_TIME("runSolver");                                                         // time keeping diagnostic

   int status = 0, mg_type = 0;                                                  // ???
   ParmParse ppMain("main");

   ppMain.query("verbosity", s_verbosity);                                       // Noisiness control ???

   // set up grids&
   Vector<DisjointBoxLayout> grids;                                              // define non-temporary structures ... ???
   Vector<ProblemDomain> domain;
   Vector<int> ref_ratio;
   Vector<Real> dx;
   int finestLevel;

   setupGrids(amrGrids, amrDomains, refRatios, amrDx, finestLevel);              // ... to create the domain and AMR blocks/boxes ???

   // initialize solver
   AMRMultiGrid<LevelData<FArrayBox> > *amrSolver;                               // Not a_amrSolver? Where is amrSolver defined?
   if ( mg_type==0 )
     {
       amrSolver = new AMRMultiGrid<LevelData<FArrayBox> >();                    // ???
     }
   else
     {
       MayDay::Error("FAS not supported");
       // int type = (int)VCYCLE;
       // ParmParse ppSolver("solver");

       // AMRFASMultiGrid<LevelData<FArrayBox> > *psol;
       // psol = new AMRFASMultiGrid<LevelData<FArrayBox> >();
       // ppSolver.query("cycle_type", type);
       // psol->setCycleType( (FASMG_type)type );
       // bool avoid_norms = false;
       // ppSolver.query("avoid_norms", avoid_norms);
       // psol->setAvoidNorms(avoid_norms);
       // int numv=1;
       // ppSolver.query("num_v_cycles", numv);
       // psol->setNumVcycles(numv);
       // amrSolver = psol;
     }
   BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;                           // Still not clear on what the BiCGStabSolver does
   bottomSolver.m_verbosity = s_verbosity-2;
   setupSolver(amrSolver, bottomSolver, amrGrids, amrDomains,
               refRatios, amrDx, finestLevel);


   // allocate solution and RHS, initialize RHS
   int numLevels = amrGrids.size();                                              // Set up containers
   Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);                          // \phi container
   Vector<LevelData<FArrayBox>* > rhs(numLevels, NULL);                          // \rho container
   // this is for convenience
   Vector<LevelData<FArrayBox>* > resid(numLevels, NULL);                        // ??? residual calculation

   for (int lev=0; lev<=finestLevel; lev++)                                      // Loop over all AMR levels
     {
       const DisjointBoxLayout& levelGrids = amrGrids[lev];                      // level dummy
       phi[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);        // create space for \phi data for each level
       rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);        // create space for \rho date for each level
       resid[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);      // create space for residual for each level
     }

   setRHS(rhs, domain, ref_ratio, dx, finestLevel );                             // set the RHS, wasn't this done in the beginning already?

   // do solve
   int iterations = 1;
   ppMain.get("iterations", iterations);                                         // ppMain = ???

   for (int iiter = 0; iiter < iterations; iiter++)                              // interate over
     {
       bool zeroInitialGuess = true;
       pout() << "about to go into solve" << endl;
       amrSolver->solve(phi, rhs, finestLevel, 0, zeroInitialGuess);             // Here is where it is all put together
       pout() << "done solve" << endl;
     }

   // write results to file

/*
*  This may be useful for testing an debugging purposes, but in the end,
*  I would prefer to keep all the information in one place with the rest of
*  the conserved/primative variables of PLUTO
*
*/

   bool writePlots = true;
   ppMain.query("writePlotFiles", writePlots);

#ifdef CH_USE_HDF5

   if (writePlots)
     {
       int numLevels = finestLevel +1;
       Vector<LevelData<FArrayBox>* > plotData(numLevels, NULL);

       pout() << "Write Plots. norm=" << amrSolver->computeAMRResidual(resid,phi,rhs,finestLevel,0) << endl;

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

   delete amrSolver;

   return status;
}

/* ********************************************************** */
/*void solveSelfGravPot(Vector<LevelData<FArrayBox>* >& a_gravpot,       // Output self-gravity potential: m_gravpot
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
/*{
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
*/

/*
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
*/

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
