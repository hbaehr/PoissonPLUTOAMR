# PoissonPLUTOAMR

### Using the elliptic solver in the Chombo libraries to solve the Poisson equation in an AMR environment

A module for solving the Poisson equation in an adaptive mesh (AMR) framework. This is a part of the PLUTO code, which through the Chombo libraries has AMR capabilities.

The majority of the code is based off a working example created by the Chombo development team and modified to work with the AMR framework existing in PLUTO.

Instructions to make this work:

-The files AMRPoissonPluto.cpp and AMRPoissonPluto.H go in $PLUTO_DIR/Src/Chombo

-Add a segment of code to AMRLevelPluto.cpp in the advance() function which calls the solver during runtime at each coarse timestep

-In the template makefile ($PLUTO_DIR/Src/Templates/makefile.chombo) add:

---INCLUDE_DIRS += -I$(CHOMBO_HOME)/src/AMRElliptic *ABOVE* the AMRTimeDependent include

---add AMRElliptic to the list of LibNames

---add AMRPoissonPluto.H and AMRPoissonPluto.o to the list of headers and object files

Additionally:

In Src/struct.h: the data structure needs to allocate new memory for the addition of Phi

In Src/pluto.h (or /Src/HD/mod_defs ???): SELFGRAV needs to be added to NVARS total, either separately or through NFLX

In Src/initialize.c: The quantity Phi needs to be added to the main data struture at initialization

In Src/MHD/rhs.c: The update sweep needs to be modified to exclude PHI

-Hans Baehr
