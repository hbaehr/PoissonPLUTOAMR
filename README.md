# PoissonPLUTOAMR

## Using the elliptic solver in the Chombo libraries to solve the Poisson equation in an AMR environment

A module for solving the Poisson equation in an adaptive mesh (AMR) framework. This is a part of the PLUTO code, which through the Chombo libraries has AMR capabilities.

The majority of the code is based off a working example created by the Chombo development team and modified to work with the AMR framework existing in PLUTO. See the example at Chombo-3.2/releasedExamples/AMRPoisson/execCell/poissonSolve.cpp

At the moment it currently works with Cartesian coordinates and periodic boundary conditions as long as additional levels of mesh refinement are not needed.

## Instructions to make this work:

### Changes necessary for the Chombo part of the Code:

-The files AMRPoissonPluto.cpp and AMRPoissonPluto.H go in $PLUTO_DIR/Src/Chombo

-Add a segment of code to AMRLevelPluto.cpp in the advance() function which calls the solver during runtime at each coarse timestep

-In the template makefile ($PLUTO_DIR/Src/Templates/makefile.chombo) add:

---INCLUDE_DIRS += -I$(CHOMBO_HOME)/src/AMRElliptic *ABOVE* the AMRTimeDependent include

---add AMRElliptic to the list of LibNames

---add AMRPoissonPluto.H and AMRPoissonPluto.o to the list of headers and object files

-add variables to transfer potential down to the patch level: AMRLevelPluto -> LevelPluto -> PatchEuler

### Additional changes necessary to the base PLUTO code:

In Src/struct.h: the data structure needs to make room for Phi

In Src/pluto.h: define the SELFGRAV variable for simple on/off switch

In Src/initialize.c: Allocate array for Phi in the main data structure

In Src/Chombo/PatchEuler.cpp: need to allocate the memory for the space in the data structure for Phi

In Src/MHD/rhs.c and Src/MHD/rhs_source.c: add the potential the source terms of the energy and momentum equations

In Src/prototypes.h: modify the functions RightHandSide and RightHandSideSource to take in the d structure

In Src/Time_Stepping/update_stage.c: modify the calls of the above functions to use the d structure

## TODO:

+Make potential available to level integrator on refined levels

+Resolve memory leak

+More sophisticated boundary conditions

+Additional refinement criteria

-Hans Baehr
