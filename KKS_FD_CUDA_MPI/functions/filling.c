#include "filling.h"

/*
 * Fills a cylinder at the location specified by the cylinder object
 */
void fillCylinder(double *phi, cylinder Cylinder,
                  domainInfo simDomain, subdomainInfo subdomain)
{
    long phaseStep = subdomain.numCompCells;
    long index;
    long xG, yG, zG;     // Global coordinates
    long x, y, z;        // Local coordinates

    double sum = 0.0;

    for (z = 1; z < subdomain.sizeZ-1; z++)
    {
        for (y = 1; y < subdomain.sizeY-1; y++)
        {
            for (x = 1; x < subdomain.sizeX-1; x++)
            {
                index = z*subdomain.zStep + y*subdomain.yStep + x;

                // Computing global coordinates
                xG = subdomain.xS + x - 1;
                yG = subdomain.yS + y - 1;
                zG = subdomain.zS + z - 1;

                if (Cylinder.phase != simDomain.numPhases-1)
                {
                    // Circle (cylinder) interior/exterior check
                    if ((Cylinder.xC - xG)*(Cylinder.xC - xG)
                        + (Cylinder.yC - yG)*(Cylinder.yC - yG)
                        <= Cylinder.radius*Cylinder.radius)
                    {
                        for (long i = 0; i < simDomain.numPhases-1; i++)
                        {
                            if (i == Cylinder.phase)
                                phi[i*phaseStep + index] = 1.0;
                            else
                                phi[i*phaseStep + index] = 0.0;
                        }
                        phi[(simDomain.numPhases-1)*phaseStep + index] = 0.0;
                    }
                    else
                    {
                        if (phi[Cylinder.phase*phaseStep + index] != 1.0)
                            phi[Cylinder.phase*phaseStep + index] = 0.0;
                    }
                }
                else
                {
                    sum = 0.0;
                    for (long i = 0; i < simDomain.numPhases-1; i++)
                    {
                        sum += phi[i*phaseStep + index];
                    }
                    if (sum <= 1.0 && sum >= 0.0)
                        phi[(simDomain.numPhases-1)*phaseStep + index] = 1.0 - sum;
                    else
                    {
                        for (long i = 0; i < simDomain.numPhases-1; i++)
                        {
                            phi[i*phaseStep + index] = 0.0;
                        }
                        phi[(simDomain.numPhases-1)*phaseStep + index] = 1.0;
                    }
                }
            }
        }
    }
}

void fillSphere(double *phi, sphere Sphere,
                domainInfo simDomain, subdomainInfo subdomain)
{
    long phaseStep = subdomain.numCompCells;
    long index;
    long xG, yG, zG;    // Global coordinates
    long x, y, z;       // Local coordinates

    double sum = 0.0;

    for (z = 1; z < subdomain.sizeZ-1; z++)
    {
        for (y = 1; y < subdomain.sizeY-1; y++)
        {
            for (x = 1; x < subdomain.sizeX-1; x++)
            {
                index = z*subdomain.zStep + y*subdomain.yStep + x;

                // Computing global coordinates
                xG = subdomain.xS + x - 1;
                yG = subdomain.yS + y - 1;
                zG = subdomain.zS + z - 1;

                if (Sphere.phase != simDomain.numPhases-1)
                {
                    // Sphere interior/exterior check
                    if ((Sphere.xC - xG)*(Sphere.xC - xG)
                        + (Sphere.yC - yG)*(Sphere.yC - yG)
                        + (Sphere.zC - zG)*(Sphere.zC - zG)
                        <= Sphere.radius*Sphere.radius)
                    {
                        for (long i = 0; i < simDomain.numPhases-1; i++)
                        {
                            if (i == Sphere.phase)
                                phi[i*phaseStep + index] = 1.0;
                            else
                                phi[i*phaseStep + index] = 0.0;
                        }
                        phi[(simDomain.numPhases-1)*phaseStep + index] = 0.0;
                    }
                    else
                    {
                        if (phi[Sphere.phase*phaseStep + index] != 1.0)
                            phi[Sphere.phase*phaseStep + index] = 0.0;
                    }
                }
                else
                {
                    sum = 0.0;
                    for (long i = 0; i < simDomain.numPhases-1; i++)
                    {
                        sum += phi[i*phaseStep + index];
                    }
                    if (sum <= 1.0 && sum >= 0.0)
                        phi[(simDomain.numPhases-1)*phaseStep + index] = 1.0 - sum;
                    else
                    {
                        for (long i = 0; i < simDomain.numPhases-1; i++)
                        {
                            phi[i*phaseStep + index] = 0.0;
                        }
                        phi[(simDomain.numPhases-1)*phaseStep + index] = 1.0;
                    }
                }
            }
        }
    }
}

void fillCube(double *phi, cube Cube,
              domainInfo simDomain, subdomainInfo subdomain)
{
    long phaseStep = subdomain.numCompCells;
    long index;
    long xG, yG, zG;     // Global coordinates

    long x, y, z;

    double sum = 0.0;

    for (z = 1; z < subdomain.sizeZ-1; z++)
    {
        for (y = 1; y < subdomain.sizeY-1; y++)
        {
            for (x = 1; x < subdomain.sizeX-1; x++)
            {
                index = z*subdomain.zStep + y*subdomain.yStep + x;

                xG = subdomain.xS + x - 1;
                yG = subdomain.yS + y - 1;
                zG = subdomain.zS + z - 1;

                if (Cube.phase != simDomain.numPhases-1)
                {
                    // Circle (cylinder) interior/exterior check
                    if (Cube.xS <= xG && Cube.xE >= xG && Cube.yS <= yG && Cube.yE >= yG && Cube.zS <= zG && Cube.zE >= zG)
                    {
                        for (long i = 0; i < simDomain.numPhases-1; i++)
                        {
                            if (i == Cube.phase)
                                phi[i*phaseStep + index] = 1.0;
                            else
                                phi[i*phaseStep + index] = 0.0;
                        }
                        phi[(simDomain.numPhases-1)*phaseStep + index] = 0.0;
                    }
                    else
                    {
                        if (phi[Cube.phase*phaseStep + index] != 1.0)
                            phi[Cube.phase*phaseStep + index] = 0.0;
                    }
                }
                else
                {
                    sum = 0.0;
                    for (long i = 0; i < simDomain.numPhases-1; i++)
                    {
                        sum += phi[i*phaseStep + index];
                    }
                    if (sum <= 1.0 && sum >= 0.0)
                        phi[(simDomain.numPhases-1)*phaseStep + index] = 1.0 - sum;
                    else
                    {
                        for (long i = 0; i < simDomain.numPhases-1; i++)
                        {
                            phi[i*phaseStep + index] = 0.0;
                        }
                        phi[(simDomain.numPhases-1)*phaseStep + index] = 1.0;
                    }
                }
            }
        }
    }
}

void fillEllipse(double *phi, ellipse Ellipse,
                 domainInfo simDomain, subdomainInfo subdomain)
{
    long phaseStep = subdomain.numCompCells;
    long index;
    long xG, yG, zG;     // Global coordinates

    long x, y, z;

    double sum = 0.0;

    for (z = 1; z < subdomain.sizeZ-1; z++)
    {
        for (y = 1; y < subdomain.sizeY-1; y++)
        {
            for (x = 1; x < subdomain.sizeX-1; x++)
            {
                index = z*subdomain.zStep + y*subdomain.yStep + x;

                xG = subdomain.xS + x - 1;
                yG = subdomain.yS + y - 1;
                zG = subdomain.zS + z - 1;

                if (Ellipse.phase != simDomain.numPhases-1)
                {
                    // Sphere interior/exterior check
                    if ((Ellipse.xC - xG)*(Ellipse.xC - xG)
                        + (Ellipse.yC - yG)*(Ellipse.yC - yG)
                        + (Ellipse.zC - zG)*(Ellipse.zC - zG)
                        <= Ellipse.major_axis*Ellipse.major_axis)
                    {
                        for (long i = 0; i < simDomain.numPhases-1; i++)
                        {
                            if (i == Ellipse.phase)
                                phi[i*phaseStep + index] = 1.0;
                            else
                                phi[i*phaseStep + index] = 0.0;
                        }
                        phi[(simDomain.numPhases-1)*phaseStep + index] = 0.0;
                    }
                    else
                    {
                        if (phi[Ellipse.phase*phaseStep + index] != 1.0)
                            phi[Ellipse.phase*phaseStep + index] = 0.0;
                    }
                }
                else
                {
                    sum = 0.0;
                    for (long i = 0; i < simDomain.numPhases-1; i++)
                    {
                        sum += phi[i*phaseStep + index];
                    }
                    if (sum <= 1.0 && sum >= 0.0)
                        phi[(simDomain.numPhases-1)*phaseStep + index] = 1.0 - sum;
                    else
                    {
                        for (long i = 0; i < simDomain.numPhases-1; i++)
                        {
                            phi[i*phaseStep + index] = 0.0;
                        }
                        phi[(simDomain.numPhases-1)*phaseStep + index] = 1.0;
                    }
                }
            }
        }
    }
}

void fillComposition(double *phi, double *comp,
                     domainInfo simDomain, subdomainInfo subdomain,
                     double ***ceq, double ***cfill)
{
    long index;
    long PHASE_FILLED = 0;
    long step = subdomain.numCompCells;
    long x, y, z;

    for (z = 1; z < subdomain.sizeZ-1; z++)
    {
        for (y = 1; y < subdomain.sizeY-1; y++)
        {
            for (x = 1; x < subdomain.sizeX-1; x++)
            {
                index = z*subdomain.zStep + y*subdomain.yStep + x;

                PHASE_FILLED = 0;

                for (long a = 0; a < simDomain.numPhases-1; a++)
                {
                    if (phi[a*step+index] == 1.0)
                    {
                        for (long b = 0; b < simDomain.numComponents-1; b++)
                        {
                            comp[b*step + index] = ceq[a][a][b];
                        }

                        PHASE_FILLED = 1;
                        break;
                    }
                }

                if (!PHASE_FILLED)
                {
                    for (long b = 0; b < simDomain.numComponents-1; b++)
                        comp[b*step + index] = cfill[simDomain.numPhases-1][simDomain.numPhases-1][b];
                }
            }
        }
    }
}

void fillDomain(domainInfo simDomain, subdomainInfo subdomain,
                simParameters simParams, double *phi, double *comp,
                fillParameters *fill)
{
    // Creating pointers for every filling type
    // This enables dynamic creation of filling objects
    sphere *Sphere;
    cylinder *Cylinder;
    cube *Cube;
    ellipse *Ellipse;

    // Maximum number of random filling-generation attempts
    long numTrials = 1e5;

    for (long i = 0; i < fill->countFill; i++)
    {
        if (fill->fillType[i] == FILLCYLINDER)
        {
            Cylinder = (cylinder*)malloc(sizeof(cylinder));

            Cylinder->xC = fill->xC[i];
            Cylinder->yC = fill->yC[i];
            Cylinder->radius = fill->radius[i];
            Cylinder->phase = fill->phase[i];
            Cylinder->zS = fill->zS[i];
            Cylinder->zE = fill->zE[i];

            if (Cylinder->phase < simDomain.numPhases)
            {
                fillCylinder(phi, *Cylinder, simDomain, subdomain);
                Cylinder->phase = simDomain.numPhases-1;
                fillCylinder(phi, *Cylinder, simDomain, subdomain);
            }

            free(Cylinder);
        }

        else if (fill->fillType[i] == FILLSPHERE)
        {
            Sphere = (sphere*)malloc(sizeof(sphere));

            Sphere->xC = fill->xC[i];
            Sphere->yC = fill->yC[i];
            Sphere->zC = fill->zC[i];
            Sphere->radius = fill->radius[i];
            Sphere->phase = fill->phase[i];

            if (Sphere->phase < simDomain.numPhases)
            {
                fillSphere(phi, *Sphere, simDomain, subdomain);
                Sphere->phase = simDomain.numPhases-1;
                fillSphere(phi, *Sphere, simDomain, subdomain);
            }
            free(Sphere);
        }

        else if (fill->fillType[i] == FILLCUBE)
        {
            Cube = (cube*)malloc(sizeof(cube));

            Cube->xS = fill->xS[i];
            Cube->xE = fill->xE[i];
            Cube->yS = fill->yS[i];
            Cube->yE = fill->yE[i];
            Cube->zS = fill->zS[i];
            Cube->zE = fill->zE[i];
            Cube->phase = fill->phase[i];

            if (Cube->phase < simDomain.numPhases)
            {
                fillCube(phi, *Cube, simDomain, subdomain);
                Cube->phase = simDomain.numPhases-1;
                fillCube(phi, *Cube, simDomain, subdomain);
            }
            free(Cube);
        }

        else if (fill->fillType[i] == FILLCYLINDERRANDOM)
        {
            double volParticle = (double)(fill->radius[i]*fill->radius[i])*M_PI;
            double volume      = (double)(simDomain.MESH_X*simDomain.MESH_Y);
            long numParticles   =  ceil(volume*fill->volFrac[i] / volParticle);

            srand48(simParams.SEED);

            Cylinder = (cylinder*)malloc(sizeof(cylinder) * numParticles);

            long count = 0;

            for (long k = 0; k < numTrials; k++)
            {
                Cylinder[count].radius = (double)fill->radius[i]*(1.0 + (drand48() - 0.5)*fill->radVar[i]);

                long distP, distM;

                // Get x-center at an adequate distance from the domain edge
                // This is done to prevent particles from coming too close as a result of the periodic B.C.
                do
                {
                    Cylinder[count].xC = simDomain.MESH_X*drand48();
                    distP = Cylinder[count].xC + Cylinder[count].radius;
                    distM = Cylinder[count].xC - Cylinder[count].radius;
                } while (!((distP < simDomain.MESH_X-1 - (fill->shieldDist[i]/2)) && (distM > fill->shieldDist[i]/2)));

                // Get y-center
                do
                {
                    Cylinder[count].yC = simDomain.MESH_Y*drand48();
                    distP = Cylinder[count].yC + Cylinder[count].radius;
                    distM = Cylinder[count].yC - Cylinder[count].radius;
                } while (!((distP < simDomain.MESH_Y-1 - (fill->shieldDist[i]/2)) && (distM > fill->shieldDist[i]/2)));

                // Force cylinder to span the entire z-axis in the domain
                Cylinder[count].zS = 0;
                Cylinder[count].zE = simDomain.MESH_Z - 1;

                Cylinder[count].phase = fill->phase[i];

                long j = 0;
                long distance, minDist;

                // Checking overlap with previously filled particles
                while (j < count)
                {
                    distance = (Cylinder[count].xC - Cylinder[j].xC)*(Cylinder[count].xC - Cylinder[j].xC)
                    + (Cylinder[count].yC - Cylinder[j].yC)*(Cylinder[count].yC - Cylinder[j].yC);

                    minDist  = ((Cylinder[count].radius + Cylinder[j].radius)*(1.0 + 0.5*fill->shieldDist[i]))
                    * ((Cylinder[count].radius + Cylinder[j].radius)*(1.0 + 0.5*fill->shieldDist[i]));

                    if (distance < minDist)
                        break;

                    j++;
                }

                if (j == count)
                {
                    if (Cylinder[count].phase < simDomain.numPhases)
                    {
                        fillCylinder(phi, Cylinder[count], simDomain, subdomain);
                        Cylinder[count].phase = simDomain.numPhases-1;
                        fillCylinder(phi, Cylinder[count], simDomain, subdomain);
                    }
                    count++;
                }

                if (count >= numParticles)
                    break;
            }

            free(Cylinder);
        }

        else if (fill->fillType[i] = FILLSPHERERANDOM)
        {
            double volParticle = (double)(fill->radius[i]*fill->radius[i]*fill->radius[i])*4.0*M_PI/3.0;
            double volume      = (double)(simDomain.MESH_X*simDomain.MESH_Y*simDomain.MESH_Z);
            long numParticles   =  ceil(volume*fill->volFrac[i] / volParticle);

            srand48(simParams.SEED);

            Sphere = (sphere*)malloc(sizeof(sphere) * numParticles);

            long count = 0;

            for (long k = 0; k < numTrials; k++)
            {
                Sphere[count].radius = (double)fill->radius[i]*(1.0 + (drand48() - 0.5)*fill->radVar[i]/100.0);

                long distP, distM;

                // Get x-center at an adequate distance from the domain edge
                // This is done to prevent particles from coming too close as a result of the periodic B.C.
                do
                {
                    Sphere[count].xC = simDomain.MESH_X*drand48();
                    distP = Sphere[count].xC + Sphere[count].radius;
                    distM = Sphere[count].xC - Sphere[count].radius;
                } while (!((distP < simDomain.MESH_X-1 - (fill->shieldDist[i]/2)) && (distM > fill->shieldDist[i]/2)));

                // Get y-center
                do
                {
                    Sphere[count].yC = simDomain.MESH_Y*drand48();
                    distP = Sphere[count].yC + Sphere[count].radius;
                    distM = Sphere[count].yC - Sphere[count].radius;
                } while (!((distP < simDomain.MESH_Y-1 - (fill->shieldDist[i]/2)) && (distM > fill->shieldDist[i]/2)));

                // Get z-center
                do
                {
                    Sphere[count].zC = simDomain.MESH_Z*drand48();
                    distP = Sphere[count].zC + Sphere[count].radius;
                    distM = Sphere[count].zC - Sphere[count].radius;
                } while (!((distP < simDomain.MESH_Z-1 - (fill->shieldDist[i]/2)) && (distM > fill->shieldDist[i]/2)));

                Sphere[count].phase = fill->phase[i];

                long j = 0;
                long distance, minDist;

                // Checking overlap with previously filled particles
                while (j < count)
                {
                    distance = (Sphere[count].xC - Sphere[j].xC)*(Sphere[count].xC - Sphere[j].xC)
                    + (Sphere[count].yC - Sphere[j].yC)*(Sphere[count].yC - Sphere[j].yC)
                    + (Sphere[count].zC - Sphere[j].zC)*(Sphere[count].zC - Sphere[j].zC);

                    minDist  = ((fill->shieldDist[i] + 1.0)*(Sphere[count].radius + Sphere[j].radius))
                    * ((fill->shieldDist[i] + 1.0)*(Sphere[count].radius + Sphere[j].radius));

                    if (distance < minDist)
                        break;

                    j++;
                }

                if (j == count)
                {
                    if (Sphere[count].phase < simDomain.numPhases)
                    {
                        fillSphere(phi, Sphere[count], simDomain, subdomain);
                        Sphere[count].phase = simDomain.numPhases-1;
                        fillSphere(phi, Sphere[count], simDomain, subdomain);
                    }
                    count++;
                }

                if (count >= numParticles)
                    break;
            }

            free(Sphere);
        }
        else
        {
            printf("Invalid filling parameters. Please check your filling file\n");
        }
    }

    fillComposition(phi, comp, simDomain, subdomain, simParams.ceq_host, simParams.cfill_host);

    // The following were allocated in fill_domain.c
    free(fill->fillType);
    free(fill->xC);
    free(fill->yC);
    free(fill->zC);
    free(fill->xS);
    free(fill->xE);
    free(fill->yS);
    free(fill->yE);
    free(fill->zS);
    free(fill->zE);
    free(fill->radius);
    free(fill->phase);
    free(fill->major_axis);
    free(fill->eccentricity);
    free(fill->rot_angle);
    free(fill->seed);
    free(fill->volFrac);
    free(fill->shieldDist);
    free(fill->radVar);
}
