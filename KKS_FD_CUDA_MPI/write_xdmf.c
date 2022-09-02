#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <arpa/inet.h>
#include <endian.h>

#define  IS_BIG_ENDIAN     (1 == htons(1))
#define  IS_LITTLE_ENDIAN  (!IS_BIG_ENDIAN)

#include "structures.h"
#include "inputReader.c"
#include "utilityFunctions.c"

domainInfo      simDomain;              // Global mesh size and cell size
controls        simControls;            // Timestep, num. of iterations, etc.
simParameters   simParams;              // All the simulation parameters and details

void write_xdmf_xml(long t, char *argv[]);

int main(int argc, char *argv[])
{
    long t, nprocs;
    long t_start, t_end;

    readInput_MPI(&simDomain, &simControls, &simParams, 0, argv);
    t_start = atol(argv[3]);
    t_end   = atol(argv[4]);

    for (t=t_start; t<=t_end; t+=simControls.saveInterval) {
        write_xdmf_xml(t, argv);
    }

    freeVars(&simDomain, &simParams);

    return 0;
}

void write_xdmf_xml(long t, char *argv[])
{
    long a, k;

    FILE *xmf = 0;

    /*
     * Open the file and write the XML description of the mesh..
     */
    char fname_write[1000];
    // 	sprintf(fname_write, "DATA/%s_%ld_%ld.xmf", argv[2], nprocs, t);
    sprintf(fname_write, "DATA/%s_%ld.xmf", argv[2], t);
    xmf = fopen(fname_write, "w");

    char fname_read[1000];
    // 	sprintf(fname_read, "%s_%ld_%ld.h5", argv[2], nprocs, t);
    sprintf(fname_read, "%s_%ld.h5", argv[2], t);

    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
    fprintf(xmf, "     <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"%ld %ld %ld\"/>\n", simDomain.MESH_Z, simDomain.MESH_Y, simDomain.MESH_X);

    fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%ld\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", simDomain.DIMENSION);
    fprintf(xmf, "       	0 0 0\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%ld\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", simDomain.DIMENSION);
    fprintf(xmf, "       	%lf %lf %lf\n", 1.0, 1.0, 1.0);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");

    for (a=0; a<simDomain.numPhases; a++) {
        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", simDomain.phaseNames[a]);
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", simDomain.MESH_Z, simDomain.MESH_Y, simDomain.MESH_X);
        fprintf(xmf, "        %s:/%s\n", fname_read, simDomain.phaseNames[a]);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
    }

    for (k=0; k<(simDomain.numComponents-1); k++) {
        fprintf(xmf, "     <Attribute Name=\"Composition_%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", simDomain.componentNames[k]);
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", simDomain.MESH_Z, simDomain.MESH_Y, simDomain.MESH_X);
        fprintf(xmf, "        %s:/Composition_%s\n", fname_read, simDomain.componentNames[k]);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
    }
    if (simControls.FUNCTION_F == 2)
    {
        for (k=0; k<(simDomain.numComponents-1); k++) {
            fprintf(xmf, "     <Attribute Name=\"Mu_%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", simDomain.componentNames[k]);
            fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", simDomain.MESH_Z, simDomain.MESH_Y, simDomain.MESH_X);
            fprintf(xmf, "        %s:/Mu_%s\n", fname_read, simDomain.componentNames[k]);
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "     </Attribute>\n");
        }
    }

    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
    return;
}


