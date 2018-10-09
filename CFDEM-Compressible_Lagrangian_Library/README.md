# CFDEM-Compressible_Lagrangian_Library

Creates liblagrangianCFDEMcomp-PUBLIC-X.Y.so

# Background
- creates a compressible lagrangian library for the CFDEM cloud 
- the source terms from CFDEM will then be in the right units for OpenFOAM compressible transport equations
- you can find the code that is modified for compressible solvers with the command

```
grep -r "compre" *
```

# Pre-requisites
- CFDEM 3.8.0
- LIGGGHTS 3.8.0

# Instructions
- place direcotry in `$CFDEM_SRC_DIR/lagrangian`
- type `wclean` then `wmake` in the cfdemParticleComp directory

