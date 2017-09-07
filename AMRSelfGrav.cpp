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

// Other things which need to be addressed:
// 1) Boundary conditions
// 2) Initial conditions
// 3) Fill the ghost zones of m_gravpot or a_gravpot by interpolation? see LevelPluto for more info
// 4) Go through each input from PLUTO and decide whether they should be  temporary 'm_' or permanent 'a_' for now they are all 'a_'

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

void ParseValue(Real* pos,                                                       // Looks like this would normally parse a BC value settings from a file
                int* dir,                                                        // dir = direction? as in dimension 1,2,3?
                Side::LoHiSide* side,                                            // Which side of the coordinate boundary?
                Real* a_values)                                                  // ??? boundary values
{
  //  ParmParse pp;
  //Real bcVal;
  //pp.get("bc_value",bcVal);                                                    // I still need to get the bc value from somewhere
  a_values[0]=0.;
}

void ParseBC(FArrayBox& a_state,                                                 // See documentation for BCfunction. State is the value of Phi within box
             const Box& a_valid,                                                 // domain of the box (might be a subdomain?) around which the BCs are 'built'
             const ProblemDomain& a_domain,                                      // This is the coarsest domain
             Real a_dx,                                                          // without grid data is this just the coarse dx
             bool a_homogeneous)                                                 // If true, ghost values are computed for a homogenous boundary condition
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

      Box valid = a_valid;                                                       // subdomain(?) around which the BCs are constructed
      for (int i=0; i<CH_SPACEDIM; ++i)                                          // CH_SPACEDIM = 2 or 3?
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))                                           // So when not periodic (in i direction) do the following
            {
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);                // One ghost cell on either side in each direction? Solver is 1st order so it makes sense
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if (!a_domain.domainBox().contains(ghostBoxLo))                    // If there are no ghost cells on the Lo side?
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
                  //       }                                                     // ... then apply Dirichlet BCs. But from where does this come?
                      DiriBC(a_state,                                            // phi values with in box valid
                             valid,                                              // subdomain box
                             a_dx,                                               // explanatory
                             true,                                               // Heh?
                             ParseValue,                                         // ???
                             i,                                                  // ???
                             Side::Lo,                                           // Lo boundary
                             1);                                                 // ???
                  //   }
                  // else
                  //   {
                  //     MayDay::Error("bogus bc flag lo");
                  //   }
                }

              if (!a_domain.domainBox().contains(ghostBoxHi))                    // So if domainBox contains no ghost cells on the Hi end?
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
   DampBC helper                                                                 // Reminder: Dirichlet conditions mean u=c at boundary, rather than du/dx=c
***********************************************************************/         // Everything beyond here still functions but this editor is weird
static void DampDiriBC( FArrayBox&      a_state,                                 // State is the value of \phi with in the box a_valid
                        const Box&      a_valid,                                 // Box around which BCs are constructed
                        const ProblemDomain& a_domain,                           // Coarsest grid domain
                        int             a_ratio,                                 // grid refinement ratio
                        int             a_dir,                                   // ???
                        Side::LoHiSide  a_side,                                  // Which side of the domain?
                        Interval&       a_interval                               // ???
                        )
{
  Real dc = 1.0, df = dc/(Real)a_ratio;                                          // ratio = 2, but what is (Real) doing there?
  for (int kk=1 ; kk <= 2 ; kk++, df += 1., dc += 1. )                           // wait so kk goes from 1 to 2 in increments of 1? df and dc increment by 1 as well?
    {
      Real fact = df/dc;                                                         // factor(?) = df/dc

      Box region = adjCellBox( a_valid, a_dir, a_side, 1 );                      // adjCellBox = ???
      region.shift( a_dir, -kk*sign(a_side) );                                   // How does this 'shift' work?

      for (BoxIterator bit(region); bit.ok(); ++bit)                             // Doing something over the 'region' but not sure what
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
void DampBC( FArrayBox& a_state,                                                 // See above for arguments
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             int a_ratio )
{

  if (!a_domain.domainBox().contains(a_state.box()))                             // If domain box does not contain a state box
    {
      if (!GlobalBCRS::s_areBCsParsed)
        {
          ParmParse pp;
          pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);                   // Get argument from file
          pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
          GlobalBCRS::s_areBCsParsed = true;
        }

      for (int i=0; i<CH_SPACEDIM; ++i)                                          // increment over all spatial dimensions
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))                                           // Does nothing if periodic in ith direction
            {
              Box valid = a_valid;                                               // Similar to above
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
                             int             a_whichpass,                        // ??? Which part of the v-cycle one is in?
                             int             a_dir,
                             Side::LoHiSide  a_side);

/***********************************************************************
   convergeGS_BC -- method to converge G-S on Dirchlet BCs                       // G-S = Gauss-Seidel method?
***********************************************************************/
void convergeGS_BC( FArrayBox& a_state,
                    const FArrayBox& a_rhs,
                    const Box& a_valid,
                    const ProblemDomain& a_domain,
                    Real a_dx,
                    int a_whichpass )                                            // See above
{

  if (!a_domain.domainBox().contains(a_state.box()))                             // Do nothing if there is no a_state?
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
*  the constant value b. The trick will be getting the right values from the
*  a_U or m_U arrays
*
*
*/

void setRHS(Vector<LevelData<FArrayBox>* > a_U,                                  // Output array of \rho: a_U[RHO]?
            Vector<ProblemDomain>& a_domain,                                     // Grid domain
            Vector<int>& m_ref_ratio,                                            // Refinement ratios between levels
            Vector<Real>& a_dx,                                                  // dx: grid spacing
            int a_finestLevel)                                                   // *** number of most refined level = len(m_level)
{
  CH_TIME("setRHS");

  for (int lev=0; lev<=a_finestLevel; lev++)                                     // Looping over all levels
    {
      LevelData<FArrayBox>& levelRhs = *(a_U(lev));                              // Each a_rhs[lev] is an FArrayBox temporarily stored as levelRhs
      const DisjointBoxLayout& levelGrids = levelRhs.getBoxes();                 // Get box indeces for each level

      // rhs is cell-centered...
      RealVect ccOffset = 0.5*a_dx[lev]*RealVect::Unit;                          // Adjust for the cell-centered location of the data // lev in PLUTO =?

      DataIterator levelDit = levelGrids.dataIterator();                         // Write data to a level
      for (levelDit.begin(); levelDit.ok(); ++levelDit)                          // Iterate over all points on the level?
        {
          FArrayBox& thisRhs = levelRhs[levelDit];                               // Dummy thisRhs for the iteration of this level

          thisRhs.setVal(0.0);                                               // Start by setting everything to 0

          BoxIterator bit(thisRhs.box());                                    // Loop over the IntVects of a Box
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();                                            // IntVect = vector of integers?, numbers corresponding to grids cells in one direction
              RealVect loc(iv);                                              // RealVect = vector of real values?, coordinate location of each cell?
              loc *= a_dx[lev];
              loc += ccOffset;

              RealVect dist = loc - center[n];                               //defining the distance of a location from the center?

              Real val = 0.0;
              thisRhs(iv,0) += a_U(iv,RHO);                                  // This looks right, but I am confused by all the different U containers
            }
         } // end loop over grids on this level
     } // end loop over levels
 }


void setupGrids(Vector<DisjointBoxLayout>& a_grids,                              // uncertain whether I will need to setup my own grid, but
                Vector<ProblemDomain>& a_domains,                                // this will be good for completeness sake
                Vector<int>& a_ref_ratio,                                        // grid parameters are already defined in other areas of
                Vector<Real>& a_dx,                                              // PLUTO Chombo and if this setup is just to define the bounds
                int& a_finestLevel)                                              // of the box then I already have that. I guess I just need to know
{                                                                                // if I should 'generate' a new grid or 'reuse' the old structure?
   CH_TIME("setupGrids");

   a_finestLevel = 0;

   // get grid generation parameters
//   int maxLevel, maxBoxSize, blockFactor;                                        // parameters to be read from file
//   Real fillRatio;                                                               // These parameters should be already provided somewhere; for now define them
   int maxLevel = 2;
   int maxBoxSize = 10000;
   int blockFactor =8;
   Real fillRatio = 0.85;

   // note that there only need to be numLevels-1 refinement ratios
   a_ref_ratios.resize(maxLevel);                                                // This confuses me, it seems like there are ususally numLevels+1 ratios
                                                                                 // AHA! maxlevel refers to the number of levels allowed, not necesarily the number used
   Vector<int>  is_periodic_int;                                                 // periodicity of the boundaries in each direction, not sure if that is the correct input
   bool is_periodic[SpaceDim];
//   ppGrids.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
   for (int dir=0; dir<SpaceDim; dir++)
     {
       is_periodic[dir] = (is_periodic_int[dir] == 1);
     }

   IntVect numCells;                                                             // All about numCells, first with initialization and then some
   Vector<int> incells(SpaceDim);                                                // parameters, but I am confused by the last line. D_DECL6???
   ppGrids.getarr("num_cells", incells, 0, SpaceDim);
   numCells = IntVect(D_DECL6(incells[0],incells[1],incells[2],
                              incells[3],incells[4],incells[5]) );

   RealVect domainSize = RealVect::Unit;                                         // Sets the physical domain of the grid
   if (ppGrids.contains("domain_size"))
     {
       Vector<Real> insize(SpaceDim);
       ppGrids.getarr("domain_size", insize, 0, SpaceDim);
       domainSize = RealVect(D_DECL6(insize[0],insize[1],insize[2],
                               insize[3],insize[4],insize[5]) );
     }

   // resize dataholders                                                         // Why do they need to be resized? For refinement?
   int maxNumLevels = maxLevel +1;
   a_grids.resize(maxNumLevels);
   a_domain.resize(maxNumLevels);
   a_dx.resize(maxNumLevels,-1);
   a_finestLevel = 0;

   // assumes dx=dy=dz                                                           // This assumption will not hold in PLUTO, but I should check
   a_dx[0] = domainSize[0]/numCells[0];                                          // where they are defined in the main code
   a_dx[1] = domainSize[1]/numCells[1];                                          // Wait is this the dx in each direction or on each level?
   a_dx[2] = domainSize[2]/numCells[2];                                          // If the later, why stop at 3 levels?

   IntVect domLo = IntVect::Zero;                                                // Defining the IntVect at the low end (first in the cycle)
   IntVect domHi  = numCells - IntVect::Unit;                                    // Defining the last IntVect

   ProblemDomain baseDomain(domLo, domHi, is_periodic);
   a_domain[0] = baseDomain;

   // set up refined domains, etc
   for (int lev=1; lev<= maxLevel; lev++)                                        // Cycling through the levels
     {                                                                           // lev=1 is first refined level
       a_domain[lev] = a_domain[lev-1];                                          // lev=2 is second refined level
       a_domain[lev].refine(a_ref_ratios[lev-1]);                                // etc.
       a_dx[lev] = a_dx[lev-1]/a_ref ratio[lev-1];
     }

   Vector<Vector<Box> > vectBoxes(maxLevel+1);                                   // No idea what vectorBoxes are, look this up !!!

   // local scope. for base-level grid generation
   {
     CH_TIME("BaseGridCreation");
     // generate base level grids

     domainSplit(baseDomain, vectBoxes[0], maxBoxSize, blockFactor);             // No idea what these are, but there seems to be some processor
                                                                                 // allotment going on
     Vector<int> procAssign(vectBoxes[0].size(), 0);

     LoadBalance(procAssign, vectBoxes[0]);

     DisjointBoxLayout baseGrids(vectBoxes[0], procAssign, baseDomain);

     a_grids[0] = baseGrids;
   }


   if (maxLevel > 0)
     {
       bool read_grids = false;
       ppGrids.query("read_in_grids", read_grids);
       if (read_grids)
         {
           for (int ilev = 1; ilev <= maxLevel; ilev++)                          // Loop over refinement levels = ilev?
             {
               const ProblemDomain& levDomain = a_domain[ilev];

               Vector<Box>   boxes;
               char boxCountVar[100];                                            // ???
               int boxCount;                                                     // Number of boxes within each level?
               sprintf(boxCountVar, "level_%d_box_count", ilev);                 // Look for this print in the output
               ppGrids.get(boxCountVar, boxCount);                               // Gets the limiting count from the parameter file
               boxes.resize(boxCount);                                           //
               for (int ibox = 0; ibox < boxCount; ibox++)                       // For each level loop through boxes until the max is reached
                 {
                   char boxLoVar[100];                                           // Why [100]?
                   char boxHiVar[100];
                   sprintf(boxLoVar, "level_%d_box_%d_lo", ilev, ibox);          // Check the output for these statements
                   sprintf(boxHiVar, "level_%d_box_%d_hi", ilev, ibox);          //
                   Vector<int> boxLo, boxHi;
                   ppGrids.getarr(boxLoVar, boxLo, 0, SpaceDim);                 // Also look in the parameter file
                   ppGrids.getarr(boxHiVar, boxHi, 0, SpaceDim);
                   IntVect ivLo(D_DECL(boxLo[0], boxLo[1], boxLo[2]));           // Here is that D_DECL again
                   IntVect ivHi(D_DECL(boxHi[0], boxHi[1], boxHi[2]));           // Is this creating the lo and hi through which the BoxIterator cycles?
                   boxes[ibox] = Box(ivLo, ivHi);
                   if (!levDomain.contains(boxes[ibox]))
                     {
                       MayDay::Error("box outside of domain");
                     }
                 }
               //check to see if level 0 domain is covered
               if (ilev == 0)
                 {
                   IntVectSet ivDom(levDomain.domainBox());
                   for (int ibox = 0; ibox < boxes.size(); ibox++)
                     {
                       ivDom -= boxes[ibox];
                     }
                   if (!ivDom.isEmpty())
                     {
                       MayDay::Error("level 0 boxes must cover the domain");
                     }
                 }
               Vector<int>  proc(boxes.size());
               LoadBalance(proc,boxes);
               a_grids[ilev] = DisjointBoxLayout(boxes, proc, levDomain);        // What is going on here ...
               a_finestLevel++;                                                  // and here?
             }

         }
       else
         {
           // tag on grad(rhs)
           int bufferSize = 1;
           BRMeshRefine meshGen(a_domain[0],                                     // So refinement is managed by some other class
                                a_ref_ratio,                                     // check TagCells.cpp to see how it is done in PLUTO. Same?
                                fillRatio,
                                blockFactor,
                                bufferSize,
                                maxBoxSize);

           // to be used by MeshRefine...
           Vector<Vector<Box> > oldMeshes(maxLevel+1);                           // How do these vector boxes relate to the a_grid structures?
           oldMeshes[0] = vectBoxes[0];
           for (int lev=1; lev<oldMeshes.size(); lev++)
             {
               oldMeshes[lev].push_back(a_domain[lev].domainBox());
             }

           Real refineThresh;
           ppGrids.get("refine_threshold", refineThresh);

           Real threshSqr = refineThresh*refineThresh;

           bool moreLevels = true;
           while (moreLevels)
             {
               // tag based on grad(rhs)                                         // I am curious to see how the grad is handled, maybe I can make use of it?
               // first need to allocate RHS
               Vector<LevelData<FArrayBox>* > tempRHS(a_finestLevel+1, NULL);
               for (int lev=0; lev<= a_finestLevel; lev++)
                 {                                                               // Looping over the levels and creating a temporary RHS from the grid layout
                   // note that we add a ghost cell to simplify gradients        // Note: a SINGLE ghost cell -- Is that what is meant by the '1'?
                   tempRHS[lev] = new LevelData<FArrayBox>(a_grids[lev],
                                                           1, IntVect::Unit);
                 }

               setRHS(tempRHS, a_domain, a_ref_ratio, a_dx,                      // This is the function that sets up the RHS, only this time on tempRHS !!!!
                      a_finestLevel);

               Vector<IntVectSet> tags(a_finestLevel+1);

               for (int lev=0; lev<a_finestLevel+1; lev++)                       // Looping over the AMR levels
                 {
                   const DisjointBoxLayout& levelGrids = a_grids[lev];
                   const LevelData<FArrayBox>& levelRHS = *tempRHS[lev];
                   IntVectSet& levelTags = tags[lev];

                   // compute mag(gradient)
                   DataIterator dit = levelGrids.dataIterator();                 // DING DING This is the iterator over grids on the level
                   for (dit.begin(); dit.ok(); ++dit)                            // So it considers that there are multiple grids per level
                     {
                       const FArrayBox& rhsFab = levelRHS[dit];
                       // local storage foer gradient
                       FArrayBox gradFab(levelGrids[dit],1);                     // So this is temporary storage?
                       gradFab.setVal(0.0);
                       Real thisGrad;

                       BoxIterator bit(levelGrids[dit]);
                       for (bit.begin(); bit.ok(); ++bit)
                         {
                           IntVect iv=bit();
                           for (int dir=0; dir<SpaceDim; dir++)
                             {
                               // use mag(undivided gradient)                    // Here is where the gradient is actually calculated
                               IntVect hi = iv + BASISV(dir);                    // Looks pretty straight forward, maybe it can be used for the potential as well
                               IntVect lo = iv - BASISV(dir);                    // It does not necessarily have to return the result to the main array, right?
                               thisGrad = rhsFab(hi,0) - rhsFab(lo,0);
                               gradFab(iv,0) += (thisGrad*thisGrad);
                             } // end loop over directions
                         } // end loop over cells

                       //gradFab now has mag(grad*dx)^2

                       // tag where mag(gradient) > tolerance^2
                       for (bit.begin(); bit.ok(); ++bit)
                         {
                           IntVect iv = bit();
                           if (gradFab(iv,0) > threshSqr)
                             {
                               levelTags |= iv;                                  // | is the bitwise inclusive or operator
                             }
                         } // end loop over cells
                     } // end loop over grids on this level

                 } // end loop over levels


               // call meshRefine.
               for (int lev=1; lev<=a_finestLevel; lev++)
                 {
                   oldMeshes[lev] = vectBoxes[lev];
                 }

               int topLevel = a_finestLevel;
               int newFinestLevel =  meshGen.regrid(vectBoxes,                   // And mesh refinement is handled by an external call
                                                    tags,                        // Check how PLUTO does this
                                                    0,
                                                    topLevel,
                                                    oldMeshes);


               // define new grids if necessary and test to see if we're done
               if (newFinestLevel > a_finestLevel)
                 {
                   a_finestLevel = newFinestLevel;

                   // setup new grid hierarchy
                   for (int lev=1; lev<=a_finestLevel; lev++)
                     {
                       Vector<int> procAssign(vectBoxes[lev].size(),0);          // More stuff I do not understand. Look into this procAssign stuff
                       LoadBalance(procAssign, vectBoxes[lev]);
                       DisjointBoxLayout levelGrids(vectBoxes[lev],
                                                    procAssign,
                                                    a_domain[lev]);
                       a_grids[lev] = levelGrids;
                     }
                   if (s_verbosity>2) pout() << "setupGrids: "<< a_finestLevel <<") size " << a_grids[a_finestLevel].size() << endl;
                 }
               else
                 {
                   moreLevels = false;
                 }

               if (a_finestLevel == maxLevel)
                 {
                   moreLevels = false;
                 }

               // clean up before starting again
               for (int lev=0; lev<tempRHS.size(); lev++)
                 {
                   delete tempRHS[lev];                                          // Empty out temporary strage
                 }

             } // end while (moreLevels)

         }

       // fill in remaining levels with empty DisjointBoxLayouts
       for (int lev= a_finestLevel+1; lev<=maxLevel; lev++)                      // So would this be the point where non-refined upper levels are filled with nulls?
         {                                                                       // No, say 5 refinement levels are possible, but only up to lev=3 is used,
           a_grids[lev] = DisjointBoxLayout();                                   // then lev = 4,5 are filled with empty boxes so they do not cause problems
         }

     }


}

/*  Setting up the parameters of the solver
*
*
*
*
*/


 void setupSolver(AMRMultiGrid<LevelData<FArrayBox> > *a_amrSolver,              // Name of the solver
             LinearSolver<LevelData<FArrayBox> >& a_bottomSolver,                // The bottom solver is what solves the Poisson equation on the coarsest level
             const Vector<DisjointBoxLayout>& a_grids,                           // Grids for each AMR level (i.e. there are multiple grids per level)
             const Vector<ProblemDomain>& a_domain,                              // Entire domain
             const Vector<int>& a_ref_ratio,                                     // Refinement ratios between levels
             const Vector<Real>& a_dx,                                           // *** dx: not sure what this value is atm
             int a_finestLevel)                                                  // *** number of most refined level
 {
   CH_TIME("setupSolver");                                                       // Timing diagnostic

//   ParmParse ppSolver("solver");                                                 // ??? WHere does ppSolver come from? Parse parameter information from input files

   int numLevels = a_finestLevel+1;                                              // 0 index count of total levels

   AMRPoissonOpFactory opFactory;                                                // Create an instance of AMRPoissonOpFactory

   // solving poisson problem here
   Real alpha =0.0;                                                              // Constants which determine form of the Poisson equation
   Real beta = 1.0/(4*3.14159265);                                               // \beta will be 1/(4*\pi*G) in cgs or code units? For now, G=1

   opFactory.define(a_domain[0],                                                 // Define the parameters that go into each instance of opFactory
                    a_grids,                                                     // Defining these parameters is the whole point of setupGrids, but these
                    a_ref_ratio,                                                 // are already defined by the main part of PLUTO-Chombo
                    a_dx[0],                                                     // so setupGrids is probably unnecessary.
                    &ParseBC, alpha, beta);

   AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;  // ??????????

   a_amrSolver->define(a_domain[0], castFact,                                    // Define
                      &a_bottomSolver, numLevels);

   // multigrid solver parameters                                                // Parameters for solving over multiple levels ???
//   int numSmooth, numMG, maxIter;                                                // So ParmParse look into the 'inputs' file and for in the instance
//   Real eps, hang;                                                               // ppSolver, extracts all values on lines which start with 'solver.'
   int numSmooth = 4;
   int numMG = 1;
   int maxIter = 100;
   Real eps = 1.0e-9;
   Real hang = 1.0e-10;
//   ppSolver.get("num_smooth", numSmooth);                                        // and returns the values for the respective variable names
//   ppSolver.get("num_mg",     numMG);
//   ppSolver.get("max_iterations", maxIter);                                      // Since I will not be reading any of this from file I will have to define it
//   ppSolver.get("tolerance", eps);
//   ppSolver.get("hang",      hang);

   Real normThresh = 1.0e-30;                                                    // What is this threshold? Look to the a_amrSolver
   a_amrSolver->setSolverParameters(numSmooth, numSmooth, numSmooth,             // Input all the parameters
                                numMG, maxIter, eps, hang, normThresh);
   a_amrSolver->m_verbosity = s_verbosity-1;

   // optional parameters
//   ppSolver.query("num_pre", a_amrSolver->m_pre);
//   ppSolver.query("num_post", a_amrSolver->m_post);
//   ppSolver.query("num_bottom", a_amrSolver->m_bottom);
 }

 int runSolver()                                                                 // Now that everything is set up (???) calculate the potential
 {
   CH_TIME("runSolver");                                                         // time keeping diagnostic

   int status = 0, mg_type = 0;                                                  // ???
//   ParmParse ppMain("main");

//   ppMain.query("verbosity", s_verbosity);                                     // Noisiness control ???
   int s_verbosity = 4;

   // set up grids&
   Vector<DisjointBoxLayout> grids;                                              // define non-temporary structures ... ???
   Vector<ProblemDomain> domain;
   Vector<int> ref_ratio;
   Vector<Real> dx;
   int finestLevel;

   setupGrids(amrGrids, amrDomains, refRatios, amrDx, finestLevel);              // but here setupGrids is called to do what? If all of these parameters are defined ...

   // initialize solver
   AMRMultiGrid<LevelData<FArrayBox> > *amrSolver;                               // Not a_amrSolver? Where is amrSolver defined? And why a pointer?
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
   BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;                           // So a v-cycle needs a direct solver on the coarsest level
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
       const DisjointBoxLayout& levelGrids = amrGrids[lev];                      // lev = level dummy
       phi[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);        // create space for \phi data for each level
       rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);        // create space for \rho date for each level
       resid[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);      // create space for residual for each level
     }

   setRHS(rhs, domain, ref_ratio, dx, finestLevel );                             // set the RHS, wasn't this done in the beginning already?
                                                                                 // or was that just the rules and parameters and this is the execution?
   // do solve
   int iterations = 1;
//   ppMain.get("iterations", iterations);                                         // ppMain = ???

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

   bool writePlots = false;
//   ppMain.query("writePlotFiles", writePlots);

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

// Output gravpot [phi] and acceleration [-grad(phi)]: used by PLUTO/Src/HD/prim_eqn.c
