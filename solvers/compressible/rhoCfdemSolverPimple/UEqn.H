// Solve the Momentum equation

// WORK AROUND
// - fvm::Sp(fvc::ddt(voidfraction),U) --> cannot figure out how to mult by rho ...
// So create new field that is voidfraction * rho, called voidfractionRho

                // Momentum predictor
                tmp<fvVectorMatrix> tUEqn
                (
                    fvm::ddt(voidfractionRho, U)  - fvm::Sp(fvc::ddt(voidfractionRho),U) 

                  + fvm::div(phi, U) - fvm::Sp(fvc::div(phi),U)

 //                 + turbulence->divDevRhoReff(U)   //comment out, replaced by term below

                  + particleCloud.divVoidfractionTau(U, voidfraction) 
                  ==
                   - fvm::Sp(Ksl,U)
                   + fvOptions(rho, U)
                );


fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();
fvOptions.constrain(UEqn);

if (momentumPredictor)
{
    if (oCorr == nOuterCorr-1)
    {

        if (modelType=="B" || modelType=="Bfull")
            solve(UEqn == - fvc::grad(p) + Ksl*Us, mesh.solver("UFinal"));
        else
            solve(UEqn == - voidfraction*fvc::grad(p) + Ksl*Us, mesh.solver("UFinal"));

        fvOptions.correct(U);
    }
    else
    {
        if (modelType=="B" || modelType=="Bfull")
            solve(UEqn == - fvc::grad(p) + Ksl*Us);
        else
            solve(UEqn == - voidfraction*fvc::grad(p) + Ksl*Us);

        fvOptions.correct(U);
    }
}

