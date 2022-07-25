#ifndef FILLING_H_
#define FILLING_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

#include "structures.h"

/*
 * Fill the domain using the input read by readFill()
 */
void fillDomain(domainInfo simDomain, subdomainInfo subdomain,
                simParameters simParams, double *phi, double *comp,
                fillParameters *fill);

/*
 * Fill a cylinder
 */
void fillCylinder(double *phi, cylinder Cylinder,
                  domainInfo simDomain, subdomainInfo subdomain);

/*
 * Fill a sphere
 */
void fillSphere(double *phi, sphere Sphere,
                domainInfo simDomain, subdomainInfo subdomain);

/*
 * Fill composition fields using the filled phase volume fraction fields
 */
void fillComposition(double *phi, double *comp,
                     domainInfo simDomain, subdomainInfo subdomain,
                     double ****ceq, double ***cfill);

#endif
