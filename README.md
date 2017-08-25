# PoissonPLUTOAMR

### Using the elliptic solver in the Chombo libraries to solve the Poisson equation in an AMR environment

My attempt to make a module for solving the Poisson equation in an adaptive mesh (AMR) framework. This is a part of the PLUTO code, which through the Chombo libraries has AMR capabilities but has not yet incorporated self-gravity to the AMR part.

The majority of the code in AMRSelfGrav.cpp is not mine, but based off a working example created by the Chombo development team and modified to work with the AMR framework existing in PLUTO.

When it is all finished there will be instructions here on how to include this with the latest PLUTO version, including necessary adjustments to the Makefiles and config files.

-Hans Baehr
