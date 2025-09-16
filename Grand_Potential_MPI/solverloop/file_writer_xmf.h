void ms_write_xmf(char * name, long step)
{
    long a, k;
	 
    /*
     * Open the file and write the XML description of the mesh..
     */

	char fname_write[1000];
    sprintf(fname_write, "DATA/%s_%ld.xmf", name, step);

    FILE  * xmf = fopen(fname_write, "w");
	
	char fname_read[1000];
    sprintf(fname_read, "%s_%ld.h5", name, step);

    fprintf(stderr, ">> writing xmf file %s\n", fname_write);
	  
    if (DIMENSION == 2)
    {
      fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
      fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
      fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
      fprintf(xmf, " <Domain>\n");
      fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
      fprintf(xmf, "     <Topology TopologyType=\"2DCoRectMesh\" NumberOfElements=\"%ld %ld\"/>\n", MESH_X+6, MESH_Y+6);
    
      fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%ld\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", DIMENSION);
      fprintf(xmf, "       	0 0\n");
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%ld\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", DIMENSION);
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
      if ((FUNCTION_F != 5) && (!GRAIN_GROWTH)) {
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
      fprintf(xmf, "       <DataItem Dimensions=\"%ld\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", DIMENSION);
      fprintf(xmf, "       	0 0 0\n");
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%ld\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", DIMENSION);
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
      if ((FUNCTION_F != 5) && (!GRAIN_GROWTH)) {  
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
}