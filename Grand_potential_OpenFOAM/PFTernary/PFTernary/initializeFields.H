
#include "clock.H"
Random randNumber(clock::getTime());


Info<< "Initialising fields... " << endl;
double xMax = max(mesh.C().component(vector::X)).value();
double yMax = max(mesh.C().component(vector::Y)).value();
double xCenter, yCenter, randTheta, randRadius;
//int NS;

//! Initial conditions for cooling simulation
if ((swcool == 1)&&(swch == 1))
{
    theta = 0.01;
    cEigenStrain = 0*cEigenStrain;
    //NS = 1;
    
    //dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    //Info << "Probability of nucleation: " << prob.value() << endl;
    
//for(int i=0; i<numSeeds; i++){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    randTheta = randNumber.globalScalar01()*(pi.value()/2);
    Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi and theta fields in seeds" << endl;
    volScalarField gaussianSeed = (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));

    //if (prob.value() > randNumber.globalScalar01()){
    theta = theta + randTheta*gaussianSeed;
    phi = phi + gaussianSeed;
    //}
//}
}

//! Initial conditions for multigrain simulation
if ((swcool == 0)&&(swch == 1))
{
    theta = 0.01;
    cEigenStrain = 0*cEigenStrain;
    
for(int i=0; i<numSeeds; i++){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    randTheta = randNumber.globalScalar01()*(pi.value()/2);
    Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi and theta fields in seeds" << endl;
    volScalarField gaussianSeed = (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));

    theta = theta + randTheta*gaussianSeed;
    phi = phi + gaussianSeed;
}
}

//! Initial conditions for precipitate growth simulation
if (swch == 2)
{
    /*theta = 0.01;
    //NS = 1;
    
    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
if (numSeeds > 1){
for(int i=0; i<numSeeds; i++){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    randRadius = (1 + randNumber.globalScalar01())*seedRadius;
    //volScalarField gaussianSeed = (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(randRadius*randRadius));
    
    volScalarField gaussianSeed = (1-phi)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) - randRadius*randRadius));

    //if (prob.value() > randNumber.globalScalar01()){
    //theta = theta + randTheta*gaussianSeed;
    phi = phi + gaussianSeed;
    //}
}
}

if (numSeeds == 1){
xCenter = 0.5*xMax;
yCenter = 0.5*yMax;
phi = phi + (1-phi)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) - seedRadius*seedRadius));//(1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));
}
}

//! Initial conditions for single dendrite solidification simulation
if (swch == 0)
{
xCenter = 0.0;
yCenter = 0.0;
phi = phi + (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));
theta = theta_val.value();
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi:" << min(phi).value() << ' ' << max(phi).value() << endl;
Info<< "Min/max theta:" << min(theta).value() << ' ' << max(theta).value() << endl;
