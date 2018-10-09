#ifdef version50
    if (modelType=="A")
    {
        volScalarField rUsed = rAU*voidfraction*rho;
        constrainPressure(p, U, phi, rUsed);
    }
    else constrainPressure(p, U, phi, rAU);
#endif
#ifndef versionExt32
#ifndef version40
    if (modelType=="A")
    {
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            #ifdef versionv1612plus
                p.boundaryFieldRef(),
            #else
                p.boundaryField(),
            #endif
            (
                phi.boundaryField()
              - (mesh.Sf().boundaryField() & U.boundaryField())
            )/(mesh.magSf().boundaryField()*rAUf.boundaryField()*voidfractionf.boundaryField())
        );
    }else
    {
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            #ifdef versionv1612plus
                p.boundaryFieldRef(),
            #else
                p.boundaryField(),
            #endif
            (
                phi.boundaryField()
              - (mesh.Sf().boundaryField() & U.boundaryField())
            )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
        );
    }
#endif
#endif
