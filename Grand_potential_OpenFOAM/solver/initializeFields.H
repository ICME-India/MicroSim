
#include "clock.H"
Random randNumber(clock::getTime());


Info<< "Initialising fields... " << endl;
double xMax = max(mesh.C().component(vector::X)).value();
double yMax = max(mesh.C().component(vector::Y)).value();
double zMax = max(mesh.C().component(vector::Z)).value();
double xCenter, yCenter, zCenter, randTheta[3], randRadius, maxvol, curvol;
//int NS;

//! For 2D simulations
if (dimensions == 2)
{
//! Initial conditions for cooling simulation
if ((swcool == 1)&&(swch == 1))
{
    theta = 0.01*vector(0,0,1);
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
    randTheta[2] = randNumber.globalScalar01()*(pi.value()/2);
    Info<< "random theta: " << randTheta[2] << endl;

    Info<< "Filling phi and theta fields in seeds" << endl;
    volScalarField gaussianSeed = (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));

    //if (prob.value() > randNumber.globalScalar01()){
    theta = theta + randTheta[2]*gaussianSeed*vector(0,0,1);
    phi = phi + gaussianSeed;
    //}
//}
}

//! Initial conditions for multigrain simulation
if ((swcool == 0)&&(swch == 1))
{
    theta = 0.01*vector(0,0,1);
    cEigenStrain = 0*cEigenStrain;
    
for(int i=0; i<numSeeds; i++){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    randTheta[2] = randNumber.globalScalar01()*(pi.value()/2);
    Info<< "random theta: " << randTheta[2] << endl;

    Info<< "Filling phi and theta fields in seeds" << endl;
    volScalarField gaussianSeed = (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));

    theta = theta + randTheta[2]*gaussianSeed*vector(0,0,1);
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
maxvol = volume_fraction*xMax*yMax;
curvol = 0;

if (numSeeds != 1){
while(curvol <= maxvol){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    randRadius = (1 + spread*(2*randNumber.globalScalar01() - 1))*seedRadius;
    curvol = curvol + pi.value()*randRadius*randRadius;
    volScalarField gaussianSeed = (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(randRadius*randRadius));
    
    //volScalarField gaussianSeed = (1-phi)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) - randRadius*randRadius));

    //if (prob.value() > randNumber.globalScalar01()){
    //theta = theta + randTheta*gaussianSeed;
    phi = phi + gaussianSeed;
    //}
}
}

if (numSeeds == 1){
xCenter = centerX;
yCenter = centerY;
phi = phi + (1-phi)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) - seedRadius*seedRadius));//(1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));
}
}

//! Initial conditions for single dendrite solidification simulation
if (swch == 0)
{
xCenter = centerX;
yCenter = centerY;
phi = phi + (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));

theta = theta_z.value()*vector(0,0,1);
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi:" << min(phi).value() << ' ' << max(phi).value() << endl;
Info<< "Min/max theta:" << min(theta).value() << ' ' << max(theta).value() << endl;
}

//! For 3D simulations
if (dimensions == 3)
{

//! Initial conditions for cooling simulation
if ((swcool == 1)&&(swch == 1))
{
    theta = 0.01*vector(1,1,1);
    cEigenStrain = 0*cEigenStrain;
    //NS = 1;

    /*dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;

for(int i=0; i<numSeeds; i++){*/
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    zCenter = randNumber.globalScalar01()*zMax;
    Info<< "xCenter, yCenter, zCenter: " << xCenter << ", " << yCenter << ", " << zCenter << endl;
    randTheta[2] = randNumber.globalScalar01()*(pi.value()/4);
    Info<< "random thetaz: " << randTheta[2] << endl;
    randTheta[0] = randNumber.globalScalar01()*(pi.value()/4);
    randTheta[1] = randNumber.globalScalar01()*(pi.value()/4);
    Info<< "random thetax: " << randTheta[0] << endl;
    Info<< "random thetay: " << randTheta[1] << endl;

    Info<< "Filling phi and theta fields in seeds" << endl;
    volScalarField gaussianSeed = (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius*seedRadius));

            //if (prob.value() > randNumber.globalScalar01()){
     theta = theta + gaussianSeed*(randTheta[0]*vector(1,0,0) + randTheta[1]*vector(0,1,0) + randTheta[2]*vector(0,0,1));
     phi = phi + gaussianSeed;

}

//! Initial conditions for multigrain simulation
if ((swcool == 0)&&(swch == 1))
{
    theta = 0.01*vector(1,1,1);
    cEigenStrain = 0*cEigenStrain;

for(int i=0; i<numSeeds; i++){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    zCenter = randNumber.globalScalar01()*zMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter, zCenter: " << xCenter << ", " << yCenter << ", " << zCenter << endl;
    randTheta[2] = randNumber.globalScalar01()*(pi.value()/4);
    Info<< "random thetaz: " << randTheta[2] << endl;
    randTheta[0] = randNumber.globalScalar01()*(pi.value()/4);
    randTheta[1] = randNumber.globalScalar01()*(pi.value()/4);
    Info<< "random thetax: " << randTheta[0] << endl;
    Info<< "random thetay: " << randTheta[1] << endl;

    Info<< "Filling phi and theta fields in seeds" << endl;
    volScalarField gaussianSeed = (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius*seedRadius));

    theta = theta + gaussianSeed*(randTheta[0]*vector(1,0,0) + randTheta[1]*vector(0,1,0) + randTheta[2]*vector(0,0,1));
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
maxvol = volume_fraction*xMax*yMax;
curvol = 0;

if (numSeeds != 1){
while(curvol <= maxvol){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    zCenter = randNumber.globalScalar01()*zMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter, zCenter: " << xCenter << ", " << yCenter << ", " << zCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    randRadius = (1 + spread*(2*randNumber.globalScalar01() - 1))*seedRadius;
    curvol = curvol + pi.value()*randRadius*randRadius;
    volScalarField gaussianSeed = (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(randRadius*randRadius));

    //volScalarField gaussianSeed = (1-phi)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter) - randRadius*randRadius));

    //if (prob.value() > randNumber.globalScalar01()){
    //theta = theta + randTheta*gaussianSeed;
    phi = phi + gaussianSeed;
    //}
}
}

if (numSeeds == 1){
xCenter = centerX;
yCenter = centerY;
zCenter = centerZ;
phi = phi + (1-phi)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter) - seedRadius*seedRadius));//(1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));
}
}

//! Initial conditions for single dendrite solidification simulation
if (swch == 0)
{
xCenter = centerX;
yCenter = centerY;
zCenter = centerZ;
phi = phi + (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius*seedRadius));

theta = theta_x.value()*vector(1,0,0) + theta_y.value()*vector(0,1,0) + theta_z.value()*vector(0,0,1);
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi:" << min(phi).value() << ' ' << max(phi).value() << endl;
Info<< "Min/max thetaz:" << min(theta.component(2)).value() << ' ' << max(theta.component(2)).value() << endl;
Info<< "Min/max thetax:" << min(theta.component(0)).value() << ' ' << max(theta.component(0)).value() << endl;
Info<< "Min/max thetay:" << min(theta.component(1)).value() << ' ' << max(theta.component(1)).value() << endl;
}
