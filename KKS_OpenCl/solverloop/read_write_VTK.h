#include<stdlib.h>
void writetofile_struct(struct variables** gridinfo,int t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  
  //Plotting in gnuplot
  for (b=0; b < NUMPHASES; b++) {
    sprintf(name,"%s_%d.vtk",Phases[b], t);
    fp=fopen(name,"w");
    write_cells_phasefield_vtk(fp, gridinfo, b, levels_of_interest);
    fclose(fp);
  }
  if (t > 0) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      sprintf(name,"%s_%d.vtk",Components[k], t);
      fp=fopen(name,"w");
      write_cells_composition_vtk(fp, gridinfo, T, k, levels_of_interest);
      fclose(fp);
    }
  } else {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      sprintf(name,"%s_%d.dat",Components[k], t);
      fp=fopen(name,"w");
      write_cells_composition_file0(fp, gridinfo, T, k, levels_of_interest);
      fclose(fp);
    }
  }
  sprintf(name,"Phaseplot_%d.vtk",t);
  fp=fopen(name,"w");
  write_cells_phaseplot_vtk(fp, gridinfo);
  fclose(fp);
#ifndef ISOTHERMAL
  sprintf(name,"Temperature_%d.vtk", t);
  fp=fopen(name,"w");
  write_cells_temperature_vtk(fp, gridinfo, levels_of_interest);
  fclose(fp);
#endif
//   for (k=0; k < NUMCOMPONENTS-1; k++) {
//     sprintf(name,"%s_%d_top_ghost_points_%d.dat",Components[k], t,levels_of_interest);
//     fp=fopen(name,"w");
//     write_cell_top_ghost_points(fp, gridinfo, T, k, levels_of_interest);
//     fclose(fp);
//   }
//   for (k=0; k < NUMCOMPONENTS-1; k++) {
//     sprintf(name,"%s_%d_bottom_ghost_points_%d.dat",Components[k], t,levels_of_interest);
//     fp=fopen(name,"w");
//     write_cell_bottom_ghost_points(fp, gridinfo, T, k, levels_of_interest);
//     fclose(fp);
//   }
}
void rdfrmfile_struct(struct variables** gridinfo,int t) {
  double x1,y1;
  long x,y;
  long gidy;
  char name[1000];
  FILE *fp;
  long b, k;
  long levels;
  
  for (b=0; b < NUMPHASES; b++) {
    sprintf(name,"%s_%d.vtk",Phases[b], t);
    fp=fopen(name,"r");
    fill_cells_phase_field_vtk(fp, gridinfo, b);
    printf("b=%ld\n",b);
    fclose(fp);
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    sprintf(name,"%s_%d.vtk",Components[k], t);
    fp=fopen(name,"r");
    fill_cells_composition_vtk(fp, gridinfo, k);
    fclose(fp);
  }
#ifndef ISOTHERMAL
  sprintf(name,"Temperature_%d.vtk",t);
  fp=fopen(name,"r");
  fill_cells_temperature_vtk(fp, gridinfo);
  fclose(fp);
#endif
#ifdef WRITECOMPOSITION
  if (t > 0) {
    compute_chemicalpotential(gridinfo);
  }
#endif
}
#ifndef ISOTHERMAL
void write_cells_temperature (FILE *fp, struct variables **gridinfo, long levels) {
long x, y,z, gidy;
  for(x=0;x<=(numx[levels]-1);x++) {
    for (y=0; y<=(numy[levels]-1); y++) {
      for (z=0;z<=(numz[levels]-1);z++) {
	gidy = x*layer_size[levels]+z*numy[levels]+y;
	fprintf(fp, "%4.3le %4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX, y*geometry[levels].DeltaY,z*geometry[levels].DeltaZ, gridinfo[levels][gidy].temperature);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
}
#endif
void write_cells_phasefield (FILE *fp, struct variables **gridinfo, long b, long levels) {
  long x, y,z, gidy;
  for(x=0;x<=(numx[levels]-1);x++) {
    for (y=0; y<=(numy[levels]-1); y++) {
      for (z=0;z<=(numz[levels]-1);z++) {
	gidy = x*layer_size[levels]+z*numy[levels]+y;
	fprintf(fp, "%4.3le %4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX, y*geometry[levels].DeltaY,z*geometry[levels].DeltaZ, gridinfo[levels][gidy].phia[b]);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
}
void write_cells_phasefield_vtk (FILE *fp, struct variables **gridinfo, long b, long levels) {
  long x, y,z, gidy;
  double phase;
  fprintf(fp,"# vtk DataFile Version 3.0\n");
  fprintf(fp,"Phasefield\n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET STRUCTURED_POINTS\n");
  fprintf(fp,"DIMENSIONS %ld %ld %ld\n",numy[0], numz[0], numx[0]);
  fprintf(fp,"ORIGIN 0 0 0\n");
  fprintf(fp,"SPACING %le %le %le\n",deltax, deltay, deltaz);
  fprintf(fp,"POINT_DATA %ld\n",numx[0]*numy[0]*numz[0]);
  fprintf(fp,"SCALARS PHASEFIELD%ld double 1\n", b);
  fprintf(fp,"LOOKUP_TABLE default\n");
  for (x=0;x<=numx[0]-1;x++) {
   for (z=0;z<=numz[0]-1;z++) {
     for (y=0; y<=numy[0]-1; y++) {
//         gidy = x*MESH_Y +y;
        gidy = x*layer_size[0]+z*numy[0]+y;
        fprintf(fp, "%le ",gridinfo[0][gidy].phia[b]);
      }
    }
    fprintf(fp, "\n");
  }
}
void write_cells_phasefield_MATRIX (FILE *fp, struct variables **gridinfo, long b, long levels) {
  long x, y, z, gidy;
  for(x=0;x<=(numx[levels]-1);x++) {
    for (y=0; y<=(numy[levels]-1); y++) {
      for (z=0;z<=(numz[levels]-1);z++) {
	gidy = x*layer_size[levels]+z*numy[levels]+y;
	fprintf(fp, "%4.3le  ",gridinfo[levels][gidy].phia[b]);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
}
void write_cells_composition (FILE *fp, struct variables **gridinfo, double T, long k, long levels) {
  long x, y,z;
  long gidy;
  long b;
  double composition;
  for(x=0; x<=(numx[levels]-1); x++) {
    for (y=0; y<=(numy[levels]-1); y++) {
      for (z=0; z<=(numz[levels]-1); z++) {
	gidy = x*layer_size[levels] + z*numy[levels] + y;
  #ifdef WRITECOMPOSITION
	composition=0.0;
  #ifdef ISOTHERMAL
	for (b=0; b < NUMPHASES; b++) {
	  composition += c_mu(gridinfo[levels][gidy].compi, T, b, k)*hphi(gridinfo[levels][gidy].phia, b);
	}
	fprintf(fp,"%4.3le %4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX, y*geometry[levels].DeltaY,z*geometry[levels].DeltaZ, composition);
  #else
	for (b=0; b < NUMPHASES; b++) {
	  composition += c_mu(gridinfo[levels][gidy].compi, gridinfo[levels][gidy].temperature, b, k)*hphi(gridinfo[levels][gidy].phia, b);
	}
	fprintf(fp,"%4.3le %4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX, y*geometry[levels].DeltaY,z*geometry[levels].DeltaZ, composition);
  #endif
  #else
	fprintf(fp,"%4.3le %4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX, y*geometry[levels].DeltaY,z*geometry[levels].DeltaZ, gridinfo[levels][gidy].compi[k]);
  #endif
      }
      fprintf(fp,"\n");
     }
    fprintf(fp,"\n");
  }
}
void write_cells_composition_file0(FILE *fp, struct variables **gridinfo, double T, long k, long levels) {
  long x, y,z,gidy;
  long b;
  double composition;
  for(x=0; x<=(numx[levels]-1); x++) {
    for (y=0; y<=(numy[levels]-1); y++) {
      for (z=0;z<=(numz[levels]-1);z++) {
	gidy = x*layer_size[levels]+z*numy[levels]+y;
	fprintf(fp,"%4.3le %4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX, y*geometry[levels].DeltaY,z*geometry[levels].DeltaZ, gridinfo[levels][gidy].compi[k]);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
}
void write_cells_composition_MATRIX(FILE *fp, struct variables **gridinfo, double T, long k, long levels) {                                                                                                               
  long x, y;
  long gidy;
  long b;
  double composition;
  for(x=0; x<=(numx[levels]-1); x++) {
    for (y=0; y<=(numy[levels]-1); y++) {
      for (z=0; z<=(numz[levels]-1); z++) {
	gidy = x*layer_size[levels]+z*numy[levels]+y;
#ifdef WRITECOMPOSITION
	composition=0.0;
#ifdef ISOTHERMAL
	for (b=0; b < NUMPHASES; b++) {
	  composition += c_mu(gridinfo[levels][gidy].compi, T, b, k)*hphi(gridinfo[levels][gidy].phia, b);
	}
	fprintf(fp,"%4.3le  ",composition);
#else
	for (b=0; b < NUMPHASES; b++) {
	  composition += c_mu(gridinfo[levels][gidy].compi, gridinfo[levels][gidy].temperature, b, k)*hphi(gridinfo[levels][gidy].phia, b);
	}
#endif
#else
	fprintf(fp,"%4.3le  ",gridinfo[levels][gidy].compi[k]);
#endif
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
}
void write_cells_composition_vtk(FILE *fp, struct variables **gridinfo, double T, long k, long levels) {
  long x, y,z;
  long gidy;
  long b;
  double composition;
  fprintf(fp,"# vtk DataFile Version 3.0\n");
  fprintf(fp,"Composition_%ld\n",k);
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET STRUCTURED_POINTS\n");
  fprintf(fp,"DIMENSIONS %ld %ld %ld\n",numy[0], numz[0], numx[0]);
  fprintf(fp,"ORIGIN 0 0 0\n");
  fprintf(fp,"SPACING %le %le %le\n",deltax, deltay, deltaz);
  fprintf(fp,"POINT_DATA %ld\n",numx[0]*numy[0]*numz[0]);
  fprintf(fp,"SCALARS COMPOSITION_%ld double 1\n",k);
  fprintf(fp,"LOOKUP_TABLE default\n");
  for (x=0;x<=numx[0]-1;x++) {
   for (z=0;z<=numz[0]-1;z++) {
     for (y=0; y<=numy[0]-1; y++) {
	gidy = x*layer_size[levels] + z*numy[levels] + y;
  #ifdef WRITECOMPOSITION
	composition=0.0;
  #ifdef ISOTHERMAL
	for (b=0; b < NUMPHASES; b++) {
	  composition += c_mu(gridinfo[levels][gidy].compi, T, b, k)*hphi(gridinfo[levels][gidy].phia, b);
	}
// 	fprintf(fp,"%4.3le %4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX, y*geometry[levels].DeltaY,z*geometry[levels].DeltaZ, composition);
        fprintf(fp, "%le ",composition);
  #else
	for (b=0; b < NUMPHASES; b++) {
	  composition += c_mu(gridinfo[levels][gidy].compi, gridinfo[levels][gidy].temperature, b, k)*hphi(gridinfo[levels][gidy].phia, b);
	}
// 	fprintf(fp,"%4.3le %4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX, y*geometry[levels].DeltaY,z*geometry[levels].DeltaZ, composition);
        fprintf(fp, "%le ",composition);
  #endif
  #else
// 	fprintf(fp,"%4.3le %4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX, y*geometry[levels].DeltaY,z*geometry[levels].DeltaZ, gridinfo[levels][gidy].compi[k]);
        fprintf(fp, "%le ",composition);
  #endif
      }
     }
    fprintf(fp,"\n");
  }
}
void write_cells_phaseplot_vtk(FILE *fp, struct variables **gridinfo) {
  long x, y, z, gidy, b;
  double phase;
  fprintf(fp,"# vtk DataFile Version 3.0\n");
  fprintf(fp,"Phaseplot\n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET STRUCTURED_POINTS\n");
  fprintf(fp,"DIMENSIONS %ld %ld %ld\n",numy[0], numz[0], numx[0]);
  fprintf(fp,"ORIGIN 0 0 0\n");
  fprintf(fp,"SPACING %le %le %le\n",deltax, deltay, deltaz);
  fprintf(fp,"POINT_DATA %ld\n",numx[0]*numy[0]*numz[0]);
  fprintf(fp,"SCALARS PHASEPLOT double 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for (x=0;x<=numx[0]-1;x++) {
   for (z=0;z<=numz[0]-1;z++) {
     for (y=0; y<=numy[0]-1; y++) {
//         gidy = x*MESH_Y +y;
        gidy = x*layer_size[0]+z*numy[0]+y;
        phase=0.0;
        for (b=0; b < NUMPHASES-1; b++) {
          phase += (b+1)*gridinfo[0][gidy].phia[b]; 
        }
        fprintf(fp, "%le ",phase);
      }
    }
    fprintf(fp, "\n");
  }
}
#ifndef ISOTHERMAL
void write_cells_temperature_vtk(FILE *fp, struct variables **gridinfo, long levels) {
  long x, y, z, gidy, b;
  double phase;
  fprintf(fp,"# vtk DataFile Version 3.0\n");
  fprintf(fp,"Temperature\n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET STRUCTURED_POINTS\n");
  fprintf(fp,"DIMENSIONS %ld %ld %ld\n",numy[0], numz[0], numx[0]);
  fprintf(fp,"ORIGIN 0 0 0\n");
  fprintf(fp,"SPACING %le %le %le\n",deltax, deltay, deltaz);
  fprintf(fp,"POINT_DATA %ld\n",numx[0]*numy[0]*numz[0]);
  fprintf(fp,"SCALARS TEMPERATURE double 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for (x=0;x<=numx[0]-1;x++) {
   for (z=0;z<=numz[0]-1;z++) {
     for (y=0; y<=numy[0]-1; y++) {
//         gidy = x*MESH_Y +y;
        gidy = x*layer_size[0]+z*numy[0]+y;
//         phase=0.0;
//         for (b=1; b < NUMPHASES-1; b++) {
//           phase += b*gridinfo[0][gidy].phia[b]; 
//         }
        fprintf(fp, "%le ",gridinfo[0][gidy].temperature);
      }
    }
    fprintf(fp, "\n");
  }
}
#endif
void write_cell_top_ghost_points(FILE *fp, struct variables **gridinfo, double T, long k, long levels) {
  long x, y,z;
  long gidy;
  long b;
  double composition;
  for(x=0; x<=(numx[levels]-1); x++) {
    for (z=0; z<=(numz[levels]-1); z++) {
      gidy = x*layer_size[levels] + z*numy[levels] + numy[levels]-1;
#ifdef WRITECOMPOSITION
      composition=0.0;
#ifdef ISOTHERMAL
      for (b=0; b < NUMPHASES; b++) {
	composition += c_mu(gridinfo[levels][gidy].compi, T, b, k)*hphi(gridinfo[levels][gidy].phia, b);
      }
      fprintf(fp,"%4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX,z*geometry[levels].DeltaZ, composition);
#else
      for (b=0; b < NUMPHASES; b++) {
	composition += c_mu(gridinfo[levels][gidy].compi, gridinfo[levels][gidy].temperature, b, k)*hphi(gridinfo[levels][gidy].phia, b);
      }
      fprintf(fp,"%4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX, z*geometry[levels].DeltaZ, composition);
#endif
#else
      fprintf(fp,"%4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX, z*geometry[levels].DeltaZ, gridinfo[levels][gidy].compi[k]);
#endif
    }
    fprintf(fp,"\n");
  }
}
void write_cell_bottom_ghost_points(FILE *fp, struct variables **gridinfo, double T, long k, long levels) {
    long x, y,z;
  long gidy;
  long b;
  double composition;
  for(x=0; x<=(numx[levels]-1); x++) {
    for (z=0; z<=(numz[levels]-1); z++) {
      gidy = x*layer_size[levels] + z*numy[levels];
#ifdef WRITECOMPOSITION
      composition=0.0;
#ifdef ISOTHERMAL
      for (b=0; b < NUMPHASES; b++) {
	composition += c_mu(gridinfo[levels][gidy].compi, T, b, k)*hphi(gridinfo[levels][gidy].phia, b);
      }
      fprintf(fp,"%4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX,z*geometry[levels].DeltaZ, composition);
#else
      for (b=0; b < NUMPHASES; b++) {
	composition += c_mu(gridinfo[levels][gidy].compi, gridinfo[levels][gidy].temperature, b, k)*hphi(gridinfo[levels][gidy].phia, b);
      }
      fprintf(fp,"%4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX,z*geometry[levels].DeltaZ, composition);
#endif
#else
      fprintf(fp,"%4.3le %4.3le %4.3le\n",(x)*geometry[levels].DeltaX,z*geometry[levels].DeltaZ, gridinfo[levels][gidy].compi[k]);
#endif
    }
    fprintf(fp,"\n");
  }
}
