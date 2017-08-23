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

/*
* A few notes:
* m_ indicates a temporary storage structure and should be used within routines
*    and re-used each timestep
*
* a_ indicates the main which stores data between timesteps
*
*
*
*/

int s_verbosity = 1;

enum probTypes {zeroRHS = 0,
                unityRHS,
                sinusoidal,
                gaussians,
                numProbTypes};

//int s_probtype = zeroRHS;
//int s_probtype = sinusoidal;
int s_probtype = gaussians;                                                      // defines the distribution of the RHS

//  -----------------------------------------
// boundary condition stuff                                                      // Are these boundary conditions between AMR levels or for the entire domain?
//  -----------------------------------------
///
/**
 */
class GlobalBCRS
{
public:
  static std::vector<bool> s_printedThatLo, s_printedThatHi;                     // ???
  static std::vector<int> s_bcLo, s_bcHi;                                        // Lo and Hi refer to what exactly?
  static RealVect s_trigvec;
  static bool s_areBCsParsed, s_valueParsed, s_trigParsed;
};

std::vector<bool> GlobalBCRS::s_printedThatLo = std::vector<bool>(SpaceDim, false); // vector<bool> is ???
std::vector<bool> GlobalBCRS::s_printedThatHi = std::vector<bool>(SpaceDim, false);
std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
RealVect          GlobalBCRS::s_trigvec = RealVect::Zero;
bool              GlobalBCRS::s_areBCsParsed= false;
bool              GlobalBCRS::s_valueParsed= false;
bool              GlobalBCRS::s_trigParsed= false;

void ParseValue(Real* pos,                                                       // Looks like this would normally parse a BC value from a file
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  //  ParmParse pp;
  //Real bcVal;
  //pp.get("bc_value",bcVal);
  a_values[0]=0.;
}

void ParseBC(FArrayBox& a_state,                                                 // This might be where I have to make some notable changes: what is a_state?
             const Box& a_valid,                                                 // a_valid is ...?
             const ProblemDomain& a_domain,                                      // This is the global domain
             Real a_dx,                                                          // without grid data is this just the coarse dx
             bool a_homogeneous)                                                 // Does this say something about the distribution? Likely FALSE, no?
{

  if (!a_domain.domainBox().contains(a_state.box()))                             // ! returns the logical negation of expression (if expr=T, !expr=F)
    {

      // if (!GlobalBCRS::s_areBCsParsed)
      //   {
      //     ParmParse pp;
      //     pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
      //     pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
      //     GlobalBCRS::s_areBCsParsed = true;
      //   }

      Box valid = a_valid;                                                       // again ??? Box valid is?
      for (int i=0; i<CH_SPACEDIM; ++i)                                          // CH_SPACEDIM = 3?
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))                                           // So when not periodic (in i direction) do the following
            {
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);                // One ghost cell on either side in each direction? Solver is 1st order so it makes sense
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if (!a_domain.domainBox().contains(ghostBoxLo))
                {
                  // if (GlobalBCRS::s_bcLo[i] == 1)
                  //   {
                  //     if (!GlobalBCRS::s_printedThatLo[i])
                  //       {
                  //         if (a_state.nComp() != 1)
                  //           {
                  //             MayDay::Error("using scalar bc function for vector");
                  //           }
                  //         GlobalBCRS::s_printedThatLo[i] = true;
                  //         if (s_verbosity>5)pout() << "const neum bcs lo for direction " << i << endl;
                  //       }
                  //     NeumBC(a_state,
                  //            valid,
                  //            a_dx,
                  //            a_homogeneous,
                  //            ParseValue,
                  //            i,
                  //            Side::Lo);
                  //   }
                  // else if (GlobalBCRS::s_bcLo[i] == 0)
                  //   {
                  //     if (!GlobalBCRS::s_printedThatLo[i])
                  //       {
                  //         if (a_state.nComp() != 1)
                  //           {
                  //             MayDay::Error("using scalar bc function for vector");
                  //           }
                  //         GlobalBCRS::s_printedThatLo[i] = true;
                  //         if (s_verbosity>5)pout() << "const diri bcs lo for direction " << i << endl;
                  //       }
                      DiriBC(a_state,                                            // Check documentation for the Dirichlet BC
                             valid,                                              // validity = consistency check?
                             a_dx,                                               // explanatory
                             true,                                               // Heh?
                             ParseValue,                                         // ???
                             i,                                                  // ???
                             Side::Lo,
                             1);                                                 // ???
                  //   }
                  // else
                  //   {
                  //     MayDay::Error("bogus bc flag lo");
                  //   }
                }

              if (!a_domain.domainBox().contains(ghostBoxHi))
                {
                  // if (GlobalBCRS::s_bcHi[i] == 1)
                  //   {
                  //     if (!GlobalBCRS::s_printedThatHi[i])
                  //       {
                  //         if (a_state.nComp() != 1)
                  //           {
                  //             MayDay::Error("using scalar bc function for vector");
                  //           }
                  //         GlobalBCRS::s_printedThatHi[i] = true;
                  //         if (s_verbosity>5)pout() << "const neum bcs hi for direction " << i << endl;
                  //       }
                  //     NeumBC(a_state,
                  //            valid,
                  //            a_dx,
                  //            a_homogeneous,
                  //            ParseValue,
                  //            i,
                  //            Side::Hi);
                  //   }
                  // else if (GlobalBCRS::s_bcHi[i] == 0)
                  //   {
                  //     if (!GlobalBCRS::s_printedThatHi[i])
                  //       {
                  //         if (a_state.nComp() != 1)
                  //           {
                  //             MayDay::Error("using scalar bc function for vector");
                  //           }
                  //         GlobalBCRS::s_printedThatHi[i] = true;
                  //         if (s_verbosity>5)pout() << "const diri bcs hi for direction " << i << endl;
                  //       }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             true,
                             ParseValue,
                             i,
                             Side::Hi,
                             1);
                  //   }
                  // else
                  //   {
                  //     MayDay::Error("bogus bc flag hi");
                  //   }
                }
            } // end if is not periodic in ith direction
        }
    }
}

//#ifdef FAS_HACKS
#if 0
/***********************************************************************
   DampBC helper
***********************************************************************/
static void DampDiriBC( FArrayBox&      a_state,
                        const Box&      a_valid,
                        const ProblemDomain& a_domain,
                        int             a_ratio,
                        int             a_dir,
                        Side::LoHiSide  a_side,
                        Interval&       a_interval
                        )
{
  Real dc = 1.0, df = dc/(Real)a_ratio;
  for (int kk=1 ; kk <= 2 ; kk++, df += 1., dc += 1. )
    {
      Real fact = df/dc;

      Box region = adjCellBox( a_valid, a_dir, a_side, 1 );
      region.shift( a_dir, -kk*sign(a_side) );

      for (BoxIterator bit(region); bit.ok(); ++bit)
        {
          const IntVect& ivTo = bit();
          for (int icomp = a_interval.begin(); icomp <= a_interval.end(); icomp++)
            {
              Real ghostVal = a_state(ivTo, icomp);
              a_state(ivTo, icomp) = fact*ghostVal;
            }
        }
    }
}
/***********************************************************************
   DampBC -- method to damp (residual) at Dirchlet BCs
***********************************************************************/
void DampBC( FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             int a_ratio )
{

  if (!a_domain.domainBox().contains(a_state.box()))
    {
      if (!GlobalBCRS::s_areBCsParsed)
        {
          ParmParse pp;
          pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
          pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
          GlobalBCRS::s_areBCsParsed = true;
        }

      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))
            {
              Box valid = a_valid;
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if ( !a_domain.domainBox().contains(ghostBoxLo) && GlobalBCRS::s_bcLo[i]==0 )
                {
                  Interval stateInterval = a_state.interval();
                  DampDiriBC(a_state,
                             valid,
                             a_domain, a_ratio,
                             i,
                             Side::Lo,
                             stateInterval
                             );
                }
              if (!a_domain.domainBox().contains(ghostBoxHi) && GlobalBCRS::s_bcHi[i]==0 )
                {
                  Interval stateInterval = a_state.interval();
                  DampDiriBC(a_state,
                             valid,
                             a_domain, a_ratio,
                             i,
                             Side::Hi,
                             stateInterval
                             );
                }
            } // end if is not periodic in ith direction
        }
    }
}

/***********************************************************************
   convergeGS_BC helper
***********************************************************************/
extern void convDiriBC_RBGS( FArrayBox&      a_state,
                             const FArrayBox& a_rhs,
                             const Box&      a_valid,
                             const ProblemDomain& a_domain,
                             Real a_dx,
                             int             a_whichpass,
                             int             a_dir,
                             Side::LoHiSide  a_side);

/***********************************************************************
   convergeGS_BC -- method to converge G-S on Dirchlet BCs
***********************************************************************/
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
          ParmParse pp;
          pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
          pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
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

/***********************************************************************
   end BC stuff
***********************************************************************/

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
