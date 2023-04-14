#include <hdf5.h>
#include<assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h> 
#include <stdbool.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "functions/global_vars.h"
#include "functions/functions.h"
#include "functions/matrix.h"
#include "functions/utility_functions.h"
#include "functions/reading_input_parameters.h"


void write_xdmf_xml(long t, char *argv[]);

void main(int argc, char *argv[]) {
  /* PASS TIME ARGUMENT ALSO: ELSE IT WILL GIVE SEGMENTATION FAULT
    */
  
  long t, nprocs;
  long t_start, t_end;
  assert(argc > 1);
//   nprocs = atol(argv[3]);
//   if (argc > 2){
//           t = atol(argv[3]);
//   }
  reading_input_parameters(argv);
  t_start = atol(argv[3]);
  t_end   = atol(argv[4]);
  
  for (t=t_start; t<=t_end; t+=saveT) {
    write_xdmf_xml(t, argv);
  }
  
  Free3M(Diffusivity, NUMPHASES, NUMCOMPONENTS-1);
  Free3M(ceq,         NUMPHASES, NUMPHASES);
  Free3M(cfill,       NUMPHASES, NUMPHASES);
  Free3M(ceq_coeffs,  NUMPHASES, NUMCOMPONENTS-1);
  Free3M(slopes,      NUMPHASES, NUMPHASES);
  Free3M(A,           NUMPHASES, NUMCOMPONENTS-1);
  Free3M(dcbdT,       NUMPHASES, NUMPHASES);

  FreeM(DELTA_T,      NUMPHASES);
  FreeM(DELTA_C,      NUMPHASES);
  FreeM(dcbdT_phase,  NUMPHASES);
  FreeM(B,            NUMPHASES);
  FreeM(Beq,          NUMPHASES);
  FreeM(dBbdT,        NUMPHASES);

  free(C);
  free(eigen_strain_phase);
  free(stiffness_phase);
  free(stiffness_t_phase);
  Free3M(cmu,NUMCOMPONENTS-1,NUMCOMPONENTS-1);
  Free3M(muc,NUMCOMPONENTS-1,NUMCOMPONENTS-1);
  Free4M(Rotation_matrix,NUMPHASES, NUMPHASES, DIMENSION);
  Free4M(Inv_Rotation_matrix,NUMPHASES, NUMPHASES, DIMENSION);
  free(Rotated_qab);
  
  return;
}
void write_xdmf_xml(long t, char *argv[]){
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
	  
    if (DIMENSION == 2) {
      fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
      fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
      fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
      fprintf(xmf, " <Domain>\n");
      fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
      fprintf(xmf, "     <Topology TopologyType=\"2DCoRectMesh\" NumberOfElements=\"%ld %ld\"/>\n", MESH_X+6, MESH_Y+6);
    
      fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", DIMENSION);
      fprintf(xmf, "       	0 0\n");
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", DIMENSION);
      fprintf(xmf, "       	%lf %lf\n", 1.0, 1.0);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");
      
      for (a=0; a<NUMPHASES; a++) {
        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Phases[a]);
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Y+6);
        fprintf(xmf, "        %s:/%s\n", fname_read,Phases[a]);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
        
      for (k=0; k<(NUMCOMPONENTS-1); k++) {
        fprintf(xmf, "     <Attribute Name=\"Mu_%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Components[k]);
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Y+6);
        fprintf(xmf, "        %s:/Mu_%s\n", fname_read,Components[k]);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }

      if (WRITECOMPOSITION) {
        for (k=0; k<(NUMCOMPONENTS-1); k++) {
          fprintf(xmf, "     <Attribute Name=\"Composition_%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Components[k]);
          fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Y+6);
          fprintf(xmf, "        %s:/Composition_%s\n", fname_read, Components[k]);
          fprintf(xmf, "       </DataItem>\n");
          fprintf(xmf, "     </Attribute>\n");
        }
      }
      if (ELASTICITY) {
        fprintf(xmf, "     <Attribute Name=\"Ux\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Y+6);
        fprintf(xmf, "        %s:/Ux\n", fname_read);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        
        fprintf(xmf, "     <Attribute Name=\"Uy\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Y+6);
        fprintf(xmf, "        %s:/Uy\n", fname_read);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      if (!ISOTHERMAL) {
        fprintf(xmf, "     <Attribute Name=\"T\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Y+6);
        fprintf(xmf, "        %s:/T\n", fname_read);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
    } else {
      fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
      fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
      fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
      fprintf(xmf, " <Domain>\n");
      fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
      fprintf(xmf, "     <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"%ld %ld %ld\"/>\n", MESH_X+6, MESH_Z+6, MESH_Y+6);
      fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", DIMENSION);
      fprintf(xmf, "       	0 0 0\n");
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", DIMENSION);
      fprintf(xmf, "       	%lf %lf %lf\n", 1.0, 1.0, 1.0);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");
      
      for (a=0; a<NUMPHASES; a++) {
        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Phases[a]);
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Z+6, MESH_Y+6);
        fprintf(xmf, "        %s:/%s\n", fname_read,Phases[a]);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
        
      for (k=0; k<(NUMCOMPONENTS-1); k++) {
        fprintf(xmf, "     <Attribute Name=\"Mu_%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Components[k]);
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Z+6, MESH_Y+6);
        fprintf(xmf, "        %s:/Mu_%s\n", fname_read,Components[k]);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }

      if (WRITECOMPOSITION) {
        for (k=0; k<(NUMCOMPONENTS-1); k++) {
          fprintf(xmf, "     <Attribute Name=\"Composition_%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Components[k]);
          fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Z+6, MESH_Y+6);
          fprintf(xmf, "        %s:/Composition_%s\n", fname_read, Components[k]);
          fprintf(xmf, "       </DataItem>\n");
          fprintf(xmf, "     </Attribute>\n");
        }
      }
      if (ELASTICITY) {
        fprintf(xmf, "     <Attribute Name=\"Ux\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Z+6, MESH_Y+6);
        fprintf(xmf, "        %s:/Ux\n", fname_read);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        
        fprintf(xmf, "     <Attribute Name=\"Uy\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Z+6, MESH_Y+6);
        fprintf(xmf, "        %s:/Uy\n", fname_read);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        
        fprintf(xmf, "     <Attribute Name=\"Uz\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Z+6, MESH_Y+6);
        fprintf(xmf, "        %s:/Uz\n", fname_read);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      if (!ISOTHERMAL) {
        fprintf(xmf, "     <Attribute Name=\"T\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", MESH_X+6, MESH_Z+6, MESH_Y+6);
        fprintf(xmf, "        %s:/T\n", fname_read);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
    }
    return;
}
// #endif // NPHASES
// #endif //DECOMP3D
// 
// 
// #endif //DIM3 

