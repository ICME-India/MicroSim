
#include "clock.H"
Random randNumber(clock::getTime());


Info<< "Initialising fields... " << endl;
double xMax = max(mesh.C().component(vector::X)).value();
double yMax = max(mesh.C().component(vector::Y)).value();
double zMax = max(mesh.C().component(vector::Z)).value();
double xCenter, yCenter, zCenter, randTheta[3], randRadius, maxvol[int(phases)-1];
int numNucl[int(phases)-1];
bool outside;
if (phases == 2)
{
//! For 2D simulations
if (dimensions == 2)
{
//! Initial conditions for cooling simulation
if ((swcool == 1)&&(swch == 1))
{
    theta = 0.01*vector(0,0,1);
    cEigenStrain1 = 0*cEigenStrain1;
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
    volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0]));

    //if (prob.value() > randNumber.globalScalar01()){
    theta += randTheta[2]*gaussianSeed*vector(0,0,1);
    phi_1 += gaussianSeed;
    //}
//}
}

//! Initial conditions for multigrain simulation
if ((swcool == 0)&&(swch == 1))
{
    theta = 0.01*vector(0,0,1);
    cEigenStrain1 = 0*cEigenStrain1;
    
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
    volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0]));

    theta += randTheta[2]*gaussianSeed*vector(0,0,1);
    phi_1 += gaussianSeed;
}
}

//! Initial conditions for precipitate growth simulation
if (swch == 2)
{
    /*theta = 0.01;
    //NS = 1;
    
    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
//curvol[0] = 0;
if (numSeeds != 1){
maxvol[0] = volumeFraction[0][0]*xMax*yMax;
numSeeds = int(maxvol[0]/(pi.value()*seedRadius[0][0]*seedRadius[0][0]));

double xcArray[int(numSeeds)], ycArray[int(numSeeds)];

xcArray[0] = randNumber.globalScalar01()*xMax;
ycArray[0] = randNumber.globalScalar01()*yMax;
Info<< "creating a new seed" << endl;
Info<< "xCenter, yCenter: " << xcArray[0] << ", " << ycArray[0] << endl;
    
for (int i=1; i<numSeeds; i++)
{
    outside = false;
    while (!outside)
    {
        xCenter = randNumber.globalScalar01()*xMax;
        yCenter = randNumber.globalScalar01()*yMax;
        for (int j=0; j<i; j++)
        {
            if (((xCenter - xcArray[j])*(xCenter - xcArray[j]) + (yCenter - ycArray[j])*(yCenter - ycArray[j])) < max(2*seedRadius+shieldDist)*max(2*seedRadius+shieldDist))
            {
                break;
            }
            if (j == i-1)
            {
                outside = true;
            }
        }
    }
    if (outside)
    {
        xcArray[i] = xCenter;
        ycArray[i] = yCenter;
        Info<< "creating a new seed" << endl;
        Info<< "xCenter, yCenter: " << xcArray[i] << ", " << ycArray[i] << endl;
    }
}


for(int i=0; i<numSeeds; i++){
//while((curvol[0]+curvol[1]+curvol[2]) <= (maxvol[0]+maxvol[1]+maxvol[2])){
    //Info<< "creating a new seed" << endl;
    //xCenter = randNumber.globalScalar01()*xMax;
    //yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    //Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    randRadius = (1 + spread[0][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[0][0];
    //if (curvol[0] <= maxvol[0]){
    //curvol[0] += pi.value()*seedRadius[0][0]*seedRadius[0][0];
    xCenter = xcArray[i];
    yCenter = ycArray[i];
    volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(randRadius*randRadius));
    
    //volScalarField gaussianSeed = (1-phi)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) - randRadius*randRadius));

    //if (prob.value() > randNumber.globalScalar01()){
    //theta = theta + randTheta*gaussianSeed;
    phi_1 += gaussianSeed;
}
}

if (numSeeds == 1){
xCenter = centerX[0][0];
yCenter = centerY[0][0];
phi_1 += (1-phi_1)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) - seedRadius[0][0]*seedRadius[0][0]));
}
}

//! Initial conditions for single dendrite solidification simulation
if (swch == 0)
{
xCenter = centerX[0][0];
yCenter = centerY[0][0];
phi_1 += (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0]));

theta = theta_z.value()*vector(0,0,1);
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;
Info<< "Min/max theta:" << min(theta).value() << ' ' << max(theta).value() << endl;

qt0 = theta.component(2);

volVectorField grad_qt0 = dimx*fvc::grad(qt0);

volScalarField grad_qt_sqr = magSqr(grad_qt0);
}

//! For 3D simulations
if (dimensions == 3)
{

//! Initial conditions for cooling simulation
if ((swcool == 1)&&(swch == 1))
{
    theta = 0.01*vector(1,1,1);
    cEigenStrain1 = 0*cEigenStrain1;
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
    volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius[0][0]*seedRadius[0][0]));

            //if (prob.value() > randNumber.globalScalar01()){
     theta += gaussianSeed*(randTheta[0]*vector(1,0,0) + randTheta[1]*vector(0,1,0) + randTheta[2]*vector(0,0,1));
     phi_1 += gaussianSeed;

}

//! Initial conditions for multigrain simulation
if ((swcool == 0)&&(swch == 1))
{
    theta = 0.01*vector(1,1,1);
    cEigenStrain1 = 0*cEigenStrain1;

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
    volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius[0][0]*seedRadius[0][0]));

    theta += gaussianSeed*(randTheta[0]*vector(1,0,0) + randTheta[1]*vector(0,1,0) + randTheta[2]*vector(0,0,1));
    phi_1 += gaussianSeed;
}
}

//! Initial conditions for precipitate growth simulation
if (swch == 2)
{
    /*theta = 0.01;
    //NS = 1;

    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
//curvol[0] = 0;

//numSeeds = 0;
//for (int i_phase = 0; i_phase < (phases-1); i_phase++){   
if (numSeeds != 1){
maxvol[0] = volumeFraction[i_phase][0]*xMax*yMax*zMax;
//curvol[i_phase] = 0;
numSeeds = int(maxvol[0]/(pi.value()*(4/3)*seedRadius[0][0]*seedRadius[0][0]*seedRadius[0][0]));
//numSeeds += numNucl[i_phase];
//}

double xcArray[int(numSeeds)], ycArray[int(numSeeds)], zcArray[int(numSeeds)];

xcArray[0] = randNumber.globalScalar01()*xMax;
ycArray[0] = randNumber.globalScalar01()*yMax;
zcArray[0] = randNumber.globalScalar01()*zMax;
Info<< "creating a new seed" << endl;
Info<< "xCenter, yCenter, zCenter: " << xcArray[0] << ", " << ycArray[0] << ", " << zcArray[0] << endl;
    
for (int i=1; i<numSeeds; i++)
{
    outside = false;
    while (!outside)
    {
        xCenter = randNumber.globalScalar01()*xMax;
        yCenter = randNumber.globalScalar01()*yMax;
        zCenter = randNumber.globalScalar01()*zMax;
        for (int j=0; j<i; j++)
        {
            if (((xCenter - xcArray[j])*(xCenter - xcArray[j]) + (yCenter - ycArray[j])*(yCenter - ycArray[j]) + (zCenter - zcArray[j])*(zCenter - zcArray[j])) < max(2*seedRadius+shieldDist)*max(2*seedRadius+shieldDist))
            {
                break;
            }
            if (j == i-1)
            {
                outside = true;
            }
        }
    }
    if (outside)
    {
        xcArray[i] = xCenter;
        ycArray[i] = yCenter;
        zcArray[i] = zCenter;
        Info<< "creating a new seed" << endl;
        Info<< "xCenter, yCenter, zCenter: " << xcArray[i] << ", " << ycArray[i] << ", " << zcArray[i] << endl;
    }
}



for(int i=0; i<numSeeds; i++){
//while(curvol[0] <= maxvol[0]){
    //Info<< "creating a new seed" << endl;
    //xCenter = randNumber.globalScalar01()*xMax;
    //yCenter = randNumber.globalScalar01()*yMax;
    //zCenter = randNumber.globalScalar01()*zMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    //Info<< "xCenter, yCenter, zCenter: " << xCenter << ", " << yCenter << ", " << zCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    randRadius = (1 + spread[0][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[0][0];
    //curvol[0] += pi.value()*(4/3)*randRadius*randRadius*randRadius;xCenter = xcArray[i];
    xCenter = xcArray[i];
    yCenter = ycArray[i];
    zCenter = zcArray[i];
    volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(randRadius*randRadius));

    //volScalarField gaussianSeed = (1-phi)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter) - randRadius*randRadius));

    //if (prob.value() > randNumber.globalScalar01()){
    //theta = theta + randTheta*gaussianSeed;
    phi_1 += gaussianSeed;
    //}
}
}

if (numSeeds == 1){
xCenter = centerX[0][0];
yCenter = centerY[0][0];
zCenter = centerZ[0][0];
phi_1 += (1-phi_1)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter) - seedRadius[0][0]*seedRadius[0][0]));//(1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));
}
}

//! Initial conditions for single dendrite solidification simulation
if (swch == 0)
{
xCenter = centerX[0][0];
yCenter = centerY[0][0];
zCenter = centerZ[0][0];
phi_1 += (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius[0][0]*seedRadius[0][0]));

theta = theta_x.value()*vector(1,0,0) + theta_y.value()*vector(0,1,0) + theta_z.value()*vector(0,0,1);
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi_1:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;
Info<< "Min/max thetaz:" << min(theta.component(2)).value() << ' ' << max(theta.component(2)).value() << endl;
Info<< "Min/max thetax:" << min(theta.component(0)).value() << ' ' << max(theta.component(0)).value() << endl;
Info<< "Min/max thetay:" << min(theta.component(1)).value() << ' ' << max(theta.component(1)).value() << endl;

qt0 = (Foam::cos(0.5*theta.component(0))*Foam::cos(0.5*theta.component(1))*Foam::cos(0.5*theta.component(2)) + Foam::sin(0.5*theta.component(0))*Foam::sin(0.5*theta.component(1))*Foam::sin(0.5*theta.component(2)));
qt1 = (Foam::sin(0.5*theta.component(0))*Foam::cos(0.5*theta.component(1))*Foam::cos(0.5*theta.component(2)) - Foam::cos(0.5*theta.component(0))*Foam::sin(0.5*theta.component(1))*Foam::sin(0.5*theta.component(2)));
qt2 = (Foam::cos(0.5*theta.component(0))*Foam::sin(0.5*theta.component(1))*Foam::cos(0.5*theta.component(2)) + Foam::sin(0.5*theta.component(0))*Foam::cos(0.5*theta.component(1))*Foam::sin(0.5*theta.component(2)));
qt3 = (Foam::cos(0.5*theta.component(0))*Foam::cos(0.5*theta.component(1))*Foam::sin(0.5*theta.component(2)) - Foam::sin(0.5*theta.component(0))*Foam::sin(0.5*theta.component(1))*Foam::cos(0.5*theta.component(2)));

volVectorField grad_qt0 = dimx*fvc::grad(qt0);
volVectorField grad_qt1 = dimx*fvc::grad(qt1);
volVectorField grad_qt2 = dimx*fvc::grad(qt2);
volVectorField grad_qt3 = dimx*fvc::grad(qt3);

volScalarField grad_qt_sqr = magSqr(grad_qt0) + magSqr(grad_qt1) + magSqr(grad_qt2) + magSqr(grad_qt3);
}
}

if (phases == 3)
{
if (dimensions == 2)
{
//! Initial conditions for precipitate growth simulation
if (swch == 2)
{
    /*theta = 0.01;
    //NS = 1;

    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
numSeeds = 0;
for (int i_phase = 0; i_phase < (phases-1); i_phase++){   
maxvol[i_phase] = volumeFraction[i_phase][0]*xMax*yMax;
//curvol[i_phase] = 0;
numNucl[i_phase] = int(maxvol[i_phase]/(pi.value()*seedRadius[i_phase][0]*seedRadius[i_phase][0]));
numSeeds += numNucl[i_phase];
}

double xcArray[int(numSeeds)], ycArray[int(numSeeds)];

xcArray[0] = randNumber.globalScalar01()*xMax;
ycArray[0] = randNumber.globalScalar01()*yMax;
Info<< "creating a new seed" << endl;
Info<< "xCenter, yCenter: " << xcArray[0] << ", " << ycArray[0] << endl;
    
for (int i=1; i<numSeeds; i++)
{
    outside = false;
    while (!outside)
    {
        xCenter = randNumber.globalScalar01()*xMax;
        yCenter = randNumber.globalScalar01()*yMax;
        for (int j=0; j<i; j++)
        {
            if (((xCenter - xcArray[j])*(xCenter - xcArray[j]) + (yCenter - ycArray[j])*(yCenter - ycArray[j])) < max(2*seedRadius+shieldDist)*max(2*seedRadius+shieldDist))
            {
                break;
            }
            if (j == i-1)
            {
                outside = true;
            }
        }
    }
    if (outside)
    {
        xcArray[i] = xCenter;
        ycArray[i] = yCenter;
        Info<< "creating a new seed" << endl;
        Info<< "xCenter, yCenter: " << xcArray[i] << ", " << ycArray[i] << endl;
    }
}

if (numSeeds != 1){
for(int i=0; i<numNucl[0]; i++){
//while((curvol[0]+curvol[1]+curvol[2]) <= (maxvol[0]+maxvol[1]+maxvol[2])){
    //Info<< "creating a new seed" << endl;
    //xCenter = randNumber.globalScalar01()*xMax;
    //yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    //Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    randRadius = (1 + spread[0][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[0][0];
    //if (curvol[0] <= maxvol[0]){
    //curvol[0] += pi.value()*seedRadius[0][0]*seedRadius[0][0];
    xCenter = xcArray[i];
    yCenter = ycArray[i];
phi_1 += (1-phi_1-phi_2)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(randRadius*randRadius));
}

//else if (curvol[1] <= maxvol[1]){
for(int i=numNucl[0]; i<numSeeds; i++){
    //curvol[1] += pi.value()*seedRadius[1][0]*seedRadius[1][0];
    randRadius = (1 + spread[1][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[1][0];
    xCenter = xcArray[i];
    yCenter = ycArray[i];
phi_2 += (1-phi_1-phi_2)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(randRadius*randRadius));
}

phi_3 = 1 - (phi_1 + phi_2);

phase_index = phi_1 + 2*phi_2 + 3*phi_3;

theta = theta_z.value()*vector(0,0,1);
}
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi_1:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;
Info<< "Min/max phi_2:" << min(phi_2).value() << ' ' << max(phi_2).value() << endl;
Info<< "Min/max phi_3:" << min(phi_3).value() << ' ' << max(phi_3).value() << endl;
Info<< "Min/max theta:" << min(theta.component(2)).value() << ' ' << max(theta.component(2)).value() << endl;
}

if (dimensions == 3)
{
//! Initial conditions for precipitate growth simulation
if (swch == 2)
{
    /*theta = 0.01;
    //NS = 1;

    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
numSeeds = 0;
for (int i_phase = 0; i_phase < (phases-1); i_phase++){   
maxvol[i_phase] = volumeFraction[i_phase][0]*xMax*yMax*zMax;
//curvol[i_phase] = 0;
numNucl[i_phase] = int(maxvol[i_phase]/(pi.value()*(4/3)*seedRadius[i_phase][0]*seedRadius[i_phase][0]*seedRadius[i_phase][0]));
numSeeds += numNucl[i_phase];
}

double xcArray[int(numSeeds)], ycArray[int(numSeeds)], zcArray[int(numSeeds)];

xcArray[0] = randNumber.globalScalar01()*xMax;
ycArray[0] = randNumber.globalScalar01()*yMax;
zcArray[0] = randNumber.globalScalar01()*zMax;
Info<< "creating a new seed" << endl;
Info<< "xCenter, yCenter, zCenter: " << xcArray[0] << ", " << ycArray[0] << ", " << zcArray[0] << endl;
    
for (int i=1; i<numSeeds; i++)
{
    outside = false;
    while (!outside)
    {
        xCenter = randNumber.globalScalar01()*xMax;
        yCenter = randNumber.globalScalar01()*yMax;
        zCenter = randNumber.globalScalar01()*zMax;
        for (int j=0; j<i; j++)
        {
            if (((xCenter - xcArray[j])*(xCenter - xcArray[j]) + (yCenter - ycArray[j])*(yCenter - ycArray[j]) + (zCenter - zcArray[j])*(zCenter - zcArray[j])) < max(2*seedRadius+shieldDist)*max(2*seedRadius+shieldDist))
            {
                break;
            }
            if (j == i-1)
            {
                outside = true;
            }
        }
    }
    if (outside)
    {
        xcArray[i] = xCenter;
        ycArray[i] = yCenter;
        zcArray[i] = zCenter;
        Info<< "creating a new seed" << endl;
        Info<< "xCenter, yCenter, zCenter: " << xcArray[i] << ", " << ycArray[i] << ", " << zcArray[i] << endl;
    }
}

if (numSeeds != 1){

for(int i=0; i<numNucl[0]; i++){
//while((curvol[0]+curvol[1]+curvol[2]) <= (maxvol[0]+maxvol[1]+maxvol[2])){
    //Info<< "creating a new seed" << endl;
    //xCenter = randNumber.globalScalar01()*xMax;
    //yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    //Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    randRadius = (1 + spread[0][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[0][0];
    //if (curvol[0] <= maxvol[0]){
    //curvol[0] += pi.value()*seedRadius[0][0]*seedRadius[0][0];
    xCenter = xcArray[i];
    yCenter = ycArray[i];
    zCenter = zcArray[i];
phi_1 += (1-phi_1-phi_2)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(randRadius*randRadius));
}

//else if (curvol[1] <= maxvol[1]){
for(int i=numNucl[0]; i<numSeeds; i++){
    //curvol[1] += pi.value()*seedRadius[1][0]*seedRadius[1][0];
    randRadius = (1 + spread[1][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[1][0];
    xCenter = xcArray[i];
    yCenter = ycArray[i];
    zCenter = zcArray[i];
phi_2 += (1-phi_1-phi_2)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(randRadius*randRadius));
}

phi_3 = 1 - (phi_1 + phi_2);

phase_index = phi_1 + 2*phi_2 + 3*phi_3;

theta = theta_x.value()*vector(1,0,0) + theta_y.value()*vector(0,1,0) + theta_z.value()*vector(0,0,1);
}
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi_1:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;
Info<< "Min/max phi_2:" << min(phi_2).value() << ' ' << max(phi_2).value() << endl;
Info<< "Min/max phi_3:" << min(phi_3).value() << ' ' << max(phi_3).value() << endl;
Info<< "Min/max thetaz:" << min(theta.component(2)).value() << ' ' << max(theta.component(2)).value() << endl;
Info<< "Min/max thetax:" << min(theta.component(0)).value() << ' ' << max(theta.component(0)).value() << endl;
Info<< "Min/max thetay:" << min(theta.component(1)).value() << ' ' << max(theta.component(1)).value() << endl;
}
}

if (phases == 4)
{
if (dimensions == 2)
{
//! Initial conditions for precipitate growth simulation
if (swch == 2)
{
    /*theta = 0.01;
    //NS = 1;

    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
numSeeds = 0;
for (int i_phase = 0; i_phase < (phases-1); i_phase++){   
maxvol[i_phase] = volumeFraction[i_phase][0]*xMax*yMax;
//curvol[i_phase] = 0;
numNucl[i_phase] = int(maxvol[i_phase]/(pi.value()*seedRadius[i_phase][0]*seedRadius[i_phase][0]));
numSeeds += numNucl[i_phase];
}

double xcArray[int(numSeeds)], ycArray[int(numSeeds)];

xcArray[0] = randNumber.globalScalar01()*xMax;
ycArray[0] = randNumber.globalScalar01()*yMax;
Info<< "creating a new seed" << endl;
Info<< "xCenter, yCenter: " << xcArray[0] << ", " << ycArray[0] << endl;
    
for (int i=1; i<numSeeds; i++)
{
    outside = false;
    while (!outside)
    {
        xCenter = randNumber.globalScalar01()*xMax;
        yCenter = randNumber.globalScalar01()*yMax;
        for (int j=0; j<i; j++)
        {
            if (((xCenter - xcArray[j])*(xCenter - xcArray[j]) + (yCenter - ycArray[j])*(yCenter - ycArray[j])) < max(2*seedRadius+shieldDist)*max(2*seedRadius+shieldDist))
            {
                break;
            }
            if (j == i-1)
            {
                outside = true;
            }
        }
    }
    if (outside)
    {
        xcArray[i] = xCenter;
        ycArray[i] = yCenter;
        Info<< "creating a new seed" << endl;
        Info<< "xCenter, yCenter: " << xcArray[i] << ", " << ycArray[i] << endl;
    }
}

if (numSeeds != 1){
for(int i=0; i<numNucl[0]; i++){
//while((curvol[0]+curvol[1]+curvol[2]) <= (maxvol[0]+maxvol[1]+maxvol[2])){
    //Info<< "creating a new seed" << endl;
    //xCenter = randNumber.globalScalar01()*xMax;
    //yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    //Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    randRadius = (1 + spread[0][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[0][0];
    //if (curvol[0] <= maxvol[0]){
    //curvol[0] += pi.value()*seedRadius[0][0]*seedRadius[0][0];
    xCenter = xcArray[i];
    yCenter = ycArray[i];
phi_1 += (1-phi_1-phi_2-phi_3)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(randRadius*randRadius));
}

//else if (curvol[1] <= maxvol[1]){
for(int i=numNucl[0]; i<numNucl[0]+numNucl[1]; i++){
    //curvol[1] += pi.value()*seedRadius[1][0]*seedRadius[1][0];
    randRadius = (1 + spread[1][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[1][0];
    xCenter = xcArray[i];
    yCenter = ycArray[i];
phi_2 += (1-phi_1-phi_2-phi_3)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(randRadius*randRadius));
}

//else if (curvol[2] <= maxvol[2]){
for(int i=numNucl[0]+numNucl[1]; i<numSeeds; i++){
    //curvol[2] += pi.value()*seedRadius[2][0]*seedRadius[2][0];
    randRadius = (1 + spread[2][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[2][0];
    xCenter = xcArray[i];
    yCenter = ycArray[i];
phi_3 += (1-phi_1-phi_2-phi_3)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(randRadius*randRadius));
}

phi_4 = 1 - (phi_1 + phi_2 + phi_3);

phase_index = phi_1 + 2*phi_2 + 3*phi_3 + 4*phi_4;

theta = theta_z.value()*vector(0,0,1);
}
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi_1:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;
Info<< "Min/max phi_2:" << min(phi_2).value() << ' ' << max(phi_2).value() << endl;
Info<< "Min/max phi_3:" << min(phi_3).value() << ' ' << max(phi_3).value() << endl;
Info<< "Min/max phi_4:" << min(phi_4).value() << ' ' << max(phi_4).value() << endl;
Info<< "Min/max theta:" << min(theta.component(2)).value() << ' ' << max(theta.component(2)).value() << endl;
}

if (dimensions == 3)
{
//! Initial conditions for precipitate growth simulation
if (swch == 2)
{
    /*theta = 0.01;
    //NS = 1;

    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
numSeeds = 0;
for (int i_phase = 0; i_phase < (phases-1); i_phase++){   
maxvol[i_phase] = volumeFraction[i_phase][0]*xMax*yMax*zMax;
//curvol[i_phase] = 0;
numNucl[i_phase] = int(maxvol[i_phase]/(pi.value()*(4/3)*seedRadius[i_phase][0]*seedRadius[i_phase][0]*seedRadius[i_phase][0]));
numSeeds += numNucl[i_phase];
}

double xcArray[int(numSeeds)], ycArray[int(numSeeds)], zcArray[int(numSeeds)];

xcArray[0] = randNumber.globalScalar01()*xMax;
ycArray[0] = randNumber.globalScalar01()*yMax;
zcArray[0] = randNumber.globalScalar01()*zMax;
Info<< "creating a new seed" << endl;
Info<< "xCenter, yCenter, zCenter: " << xcArray[0] << ", " << ycArray[0] << ", " << zcArray[0] << endl;
    
for (int i=1; i<numSeeds; i++)
{
    outside = false;
    while (!outside)
    {
        xCenter = randNumber.globalScalar01()*xMax;
        yCenter = randNumber.globalScalar01()*yMax;
        zCenter = randNumber.globalScalar01()*zMax;
        for (int j=0; j<i; j++)
        {
            if (((xCenter - xcArray[j])*(xCenter - xcArray[j]) + (yCenter - ycArray[j])*(yCenter - ycArray[j]) + (zCenter - zcArray[j])*(zCenter - zcArray[j])) < max(2*seedRadius+shieldDist)*max(2*seedRadius+shieldDist))
            {
                break;
            }
            if (j == i-1)
            {
                outside = true;
            }
        }
    }
    if (outside)
    {
        xcArray[i] = xCenter;
        ycArray[i] = yCenter;
        zcArray[i] = zCenter;
        Info<< "creating a new seed" << endl;
        Info<< "xCenter, yCenter, zCenter: " << xcArray[i] << ", " << ycArray[i] << ", " << zcArray[i] << endl;
    }
}

if (numSeeds != 1){

for(int i=0; i<numNucl[0]; i++){
//while((curvol[0]+curvol[1]+curvol[2]) <= (maxvol[0]+maxvol[1]+maxvol[2])){
    //Info<< "creating a new seed" << endl;
    //xCenter = randNumber.globalScalar01()*xMax;
    //yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    //Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    randRadius = (1 + spread[0][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[0][0];
    //if (curvol[0] <= maxvol[0]){
    //curvol[0] += pi.value()*seedRadius[0][0]*seedRadius[0][0];
    xCenter = xcArray[i];
    yCenter = ycArray[i];
    zCenter = zcArray[i];
phi_1 += (1-phi_1-phi_2-phi_3)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(randRadius*randRadius));
}

//else if (curvol[1] <= maxvol[1]){
for(int i=numNucl[0]; i<numNucl[0]+numNucl[1]; i++){
    //curvol[1] += pi.value()*seedRadius[1][0]*seedRadius[1][0];
    randRadius = (1 + spread[1][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[1][0];
    xCenter = xcArray[i];
    yCenter = ycArray[i];
    zCenter = zcArray[i];
phi_2 += (1-phi_1-phi_2-phi_3)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(randRadius*randRadius));
}

//else if (curvol[2] <= maxvol[2]){
for(int i=numNucl[0]+numNucl[1]; i<numSeeds; i++){
    //curvol[2] += pi.value()*seedRadius[2][0]*seedRadius[2][0];
    randRadius = (1 + spread[2][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[2][0];
    xCenter = xcArray[i];
    yCenter = ycArray[i];
    zCenter = zcArray[i];
phi_3 += (1-phi_1-phi_2-phi_3)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(randRadius*randRadius));
}

phi_4 = 1 - (phi_1 + phi_2 + phi_3);

phase_index = phi_1 + 2*phi_2 + 3*phi_3 + 4*phi_4;

theta = theta_x.value()*vector(1,0,0) + theta_y.value()*vector(0,1,0) + theta_z.value()*vector(0,0,1);
}
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi_1:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;
Info<< "Min/max phi_2:" << min(phi_2).value() << ' ' << max(phi_2).value() << endl;
Info<< "Min/max phi_3:" << min(phi_3).value() << ' ' << max(phi_3).value() << endl;
Info<< "Min/max phi_4:" << min(phi_4).value() << ' ' << max(phi_4).value() << endl;
Info<< "Min/max thetaz:" << min(theta.component(2)).value() << ' ' << max(theta.component(2)).value() << endl;
Info<< "Min/max thetax:" << min(theta.component(0)).value() << ' ' << max(theta.component(0)).value() << endl;
Info<< "Min/max thetay:" << min(theta.component(1)).value() << ' ' << max(theta.component(1)).value() << endl;
}
}

