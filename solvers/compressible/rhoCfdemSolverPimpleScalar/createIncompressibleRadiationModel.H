    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(Temp)
    );

    dimensionedScalar rhoCpRef
    (
        "rhoCpRef",
        dimEnergy/dimMass/dimTemperature,
        1.0
    );

    if (radiation->radiation())
    {
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false // Do not register
            )
        );

        dimensionedScalar rhoRef
        (
            "rhoRef",
            dimDensity,
            transportProperties
        );

        dimensionedScalar CpRef
        (
            "CpRef",
            dimSpecificHeatCapacity,
            transportProperties
        );

        rhoCpRef = CpRef;
    }
