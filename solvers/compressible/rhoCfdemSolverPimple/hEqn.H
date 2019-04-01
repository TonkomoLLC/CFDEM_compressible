{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(voidfractionRho, he) + fvm::div(phi, he)
      + fvc::ddt(voidfractionRho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -DpDt
        )
      - fvm::laplacian(turbulence->alphaEff()*voidfraction, he)
     ==
        voidfractionRho*(U&g)
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
}
