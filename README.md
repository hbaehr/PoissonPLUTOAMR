# PoissonPLUTOAMR

### Using the elliptic solver in the Chombo libraries to solve the Poisson equation in an AMR environment

A module for solving the Poisson equation in an adaptive mesh (AMR) framework. This is a part of the PLUTO code, which through the Chombo libraries has AMR capabilities.

The majority of the code is based off a working example created by the Chombo development team and modified to work with the AMR framework existing in PLUTO.

Instructions to make this work:

-The files AMRPoissonPluto.cpp and AMRPoissonPluto.H go in $PLUTO_DIR/Src/Chombo

-Add a segment of code to AMRLevelPluto.cpp in the advance() function which calls the solver during runtime at each coarse timestep

-Add the Self_Gravity directory to the main PLUTO/Src directory: This part is still a work in progress

-In the template makefile ($PLUTO_DIR/Src/Templates/makefile.chombo) add:

---INCLUDE_DIRS += -I$(CHOMBO_HOME)/src/AMRElliptic *ABOVE* the AMRTimeDependent include

---add AMRElliptic to the list of LibNames

---add AMRPoissonPluto.H and AMRPoissonPluto.o to the list of headers and object files

---add variables to transfer potential down to the patch level: AMRLevelPluto -> LevelPluto -> PatchEuler

Additionally:

In Src/struct.h: the data structure needs to make room for Phi

In Src/Chombo/PatchEuler.cpp: need to allocate the memory for the space in the data structure for Phi

TODO:

-Make potential available to level integrator on refined levels

-Resolve memory leak

-More sophisticated boundary conditions

-Additional refinedment criterion

-Hans Baehr
