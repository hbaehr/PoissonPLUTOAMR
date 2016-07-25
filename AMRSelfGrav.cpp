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

#include "UsingNamespace.H"

#ifdef CH_LINUX
// Should be undefined by default
#define TRAP_FPE
#undef  TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

/* List of things which need to happen
 *
 * Read in grid/mesh structure
 *  - Identify areas with certain levels of refinement
 *  - Communication across level boundaries?
 * Read in grid/mesh data
 * Calculate self-gravitational potential
 * Calculate aceleration due to self-gravity
 * Export calculations to main part of PLUTO
 *  - Store calculations somewhere (temporarily)?
 * 
 *
 */

// Read in grid/mesh structure
    
    // How to calculate for different AMR Levels? Start with calculation on most refined level
    // and then move to the next level while excluding all higher levels
    
// Read in grid/mesh data
    
    // Primitive variables from the input conservative array

    
// Use CHOMBO functions to solve for the self-gravitational potential
    
AMRPoissonOpFactory::define(a_coarseDomain   // domain at the coarsest level
                            a_grids          // AMR heirarchy
                            a_refRatios      // refinement ratios between levels
                            a_coarsedx       // grid spacing at the coarsest level
                            a_bc             // boundary conditions
                            a_alpha=1.0      // identity coefficient = 4*pi*G
                            a_beta=1.0       // Lapacian coefficient
                            );
        
// Calculate physics of self-gravity grad(phi) = g

 // Find the gradient operator in the Chombo library

// Export calculations to the main part of PLUTO
 
 // 

// Clean up and close up shop