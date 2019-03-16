    if (modelType=="A")
    {
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
                p.boundaryFieldRef(),
            (
                phi.boundaryField()
              - (mesh.Sf().boundaryField() & U.boundaryField())
            )/(mesh.magSf().boundaryField()*rAUf.boundaryField()*voidfractionf.boundaryField())
        );
    }else
    {
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
                p.boundaryFieldRef(),
            (
                phi.boundaryField()
              - (mesh.Sf().boundaryField() & U.boundaryField())
            )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
        );
    }

