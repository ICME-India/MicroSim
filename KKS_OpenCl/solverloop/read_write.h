void writetofile_struct(struct variables* gridinfo,int t) {
  long x,y;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  
  //Plotting in gnuplot
  for (b=0; b < NUMPHASES; b++) {
    if (b == (NUMPHASES-1)) {
      sprintf(name,"%s_%d.dat",Phases[b], t);
      fp=fopen(name,"w");
      write_cells_phasefield(fp, gridinfo, b);
      fclose(fp);
     }
  }
  if (t > 0) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      sprintf(name,"%s_%d.dat",Components[k], t);
      fp=fopen(name,"w");
      write_cells_composition(fp, gridinfo, T, k);
      fclose(fp);
    }
  } else {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      sprintf(name,"%s_%d.dat",Components[k], t);
      fp=fopen(name,"w");
      write_cells_composition_file0(fp, gridinfo, T, k);
      fclose(fp);
    }
  }
#ifndef ISOTHERMAL
  sprintf(name,"Temperature_%d.dat", t);
  fp=fopen(name,"w");
  write_cells_temperature(fp, gridinfo);
  fclose(fp);
#endif
}
void rdfrmfile_struct(struct variables* gridinfo,int t) {
  double x1,y1;
  long x,y;
  long gidy;
  char name[1000];
  FILE *fp;
  long b, k;
  for (b=0; b < NUMPHASES; b++) {
    if (b==NUMPHASES-1) {
     sprintf(name,"%s_%d.dat",Phases[b], t);
     fp=fopen(name,"r");
     fill_cells_phase_field_liquid(fp, gridinfo, b);
     printf("b=%ld\n",b);
     fclose(fp);
    }
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    sprintf(name,"%s_%d.dat",Components[k], t);
    fp=fopen(name,"r");
    fill_cells_composition(fp, gridinfo, k);
    fclose(fp);
  }
#ifdef WRITECOMPOSITION
  if (t > 0) {
    compute_chemicalpotential(gridinfo);
  }
#endif
}
#ifndef ISOTHERMAL
void write_cells_temperature (FILE *fp, struct variables *gridinfo) {
  long x, y, gidy;
  for(x=0;x<=(MESH_X-1);x++) {
    for (y=0; y<=(MESH_Y-1); y++) {
      gidy = x*MESH_Y +y;
      fprintf(fp, "%4.3le %4.3le %4.3le\n",(x)*deltax, y*deltay, gridinfo[gidy].temperature);
    }
    fprintf(fp,"\n");
  }
}
#endif
void write_cells_phasefield (FILE *fp, struct variables *gridinfo, long b) {
  long x, y, gidy;
  for(x=0;x<=(MESH_X-1);x++) {
    for (y=0; y<=(MESH_Y-1); y++) {
      gidy = x*MESH_Y +y;
      fprintf(fp, "%4.3le %4.3le %4.3le\n",(x)*deltax, y*deltay, gridinfo[gidy].phia[b]);
    }
    fprintf(fp,"\n");
  }
}
void write_cells_phasefield_MATRIX (FILE *fp, struct variables *gridinfo, long b) {
  long x, y, gidy;
  for(x=0;x<=(MESH_X-1);x++) {
    for (y=0; y<=(MESH_Y-1); y++) {
      gidy = x*MESH_Y +y;
      fprintf(fp, "%4.3le  ",gridinfo[gidy].phia[b]);
    }
    fprintf(fp,"\n");
  }
}
void write_cells_composition (FILE *fp, struct variables *gridinfo, double T, long k) {
  long x, y;
  long gidy;
  long b;
  double composition;
  
  for(x=0; x<=(MESH_X-1); x++) {
    for (y=0; y<=(MESH_Y-1); y++) {
      gidy = x*MESH_Y +y;
#ifdef WRITECOMPOSITION
      composition=0.0;
#ifdef ISOTHERMAL
      for (b=0; b < NUMPHASES; b++) {
        composition += c_mu(gridinfo[gidy].compi, T, b, k)*hphi(gridinfo[gidy].phia, b);
      }
      fprintf(fp,"%4.3le %4.3le %4.3le\n",(x)*deltax, y*deltay, composition);
#else
      for (b=0; b < NUMPHASES; b++) {
        composition += c_mu(gridinfo[gidy].compi, gridinfo[gidy].temperature, b, k)*hphi(gridinfo[gidy].phia, b);
      }
      fprintf(fp,"%4.3le %4.3le %4.3le\n",(x)*deltax, y*deltay, composition);
#endif
#else
      fprintf(fp,"%4.3le %4.3le %4.3le\n",(x)*deltax, y*deltay, gridinfo[gidy].compi[k]);
#endif
    }
    fprintf(fp,"\n");
  }
}
void write_cells_composition_file0(FILE *fp, struct variables *gridinfo, double T, long k) {
  long x, y;
  long gidy;
  long b;
  double composition;
  
  for(x=0; x<=(MESH_X-1); x++) {
    for (y=0; y<=(MESH_Y-1); y++) {
      gidy = x*MESH_Y +y;
      fprintf(fp,"%4.3le %4.3le %4.3le\n",(x)*deltax, y*deltay, gridinfo[gidy].compi[k]);
    }
    fprintf(fp,"\n");
  }
}
void write_cells_composition_MATRIX(FILE *fp, struct variables *gridinfo, double T, long k) {
  long x, y;
  long gidy;
  long b;
  double composition;
  
  for(x=0; x<=(MESH_X-1); x++) {
    for (y=0; y<=(MESH_Y-1); y++) {
      gidy = x*MESH_Y +y;
#ifdef WRITECOMPOSITION
      composition=0.0;
#ifdef ISOTHERMAL
      for (b=0; b < NUMPHASES; b++) {
        composition += c_mu(gridinfo[gidy].compi, T, b, k)*hphi(gridinfo[gidy].phia, b);
      }
      fprintf(fp,"%4.3le  ",composition);
#else
       for (b=0; b < NUMPHASES; b++) {
        composition += c_mu(gridinfo[gidy].compi, gridinfo[gidy].temperature, b, k)*hphi(gridinfo[gidy].phia, b);
      }
#endif
#else
      fprintf(fp,"%4.3le  ",gridinfo[gidy].compi[k]);
#endif
    }
    fprintf(fp,"\n");
  }
}
// void fill_cells_phase_field (FILE *fp, struct variables* gridinfo, long b) {
//   long x, y;
//   double x1, y1;
//   long gidy;
//   long a;
//   double sum;
//   for(x=0;x<=(MESH_X-1);x++) {
//     for (y=0; y<=(MESH_Y-1); y++) {
//       gidy = x*MESH_Y + y;
//       sum = 0.0;
//       if (b==(NUMPHASES-1)) {
//         for (a=0; a < NUMPHASES-1; a++) {
//           sum += gridinfo[gidy].phia[a];
//         }
//         gridinfo[gidy].phia[b] = 1.0 - sum;
//       } else {
//        fscanf(fp, "%le %le %le\n",&x1,&y1,&gridinfo[gidy].phia[b]);
//       }
//     }
//   }
// }
void fill_cells_phase_field_liquid(FILE *fp, struct variables* gridinfo, long b) {
  long x, y;
  double x1, y1;
  long gidy;
  long a;
  double sum;
  for(x=0;x<=(MESH_X-1);x++) {
    for (y=0; y<=(MESH_Y-1); y++) {
      gidy = x*MESH_Y + y;
      sum = 0.0;
      if (b==(NUMPHASES-1)) {
       fscanf(fp, "%le %le %le\n",&x1,&y1,&gridinfo[gidy].phia[b]);
        gridinfo[gidy].phia[0] = 1.0 - gridinfo[gidy].phia[NUMPHASES-1];
      }
    }
  }
}
void fill_cells_composition(FILE *fp, struct variables* gridinfo, long k) {
  long x, y;
  double x1, y1;
  long gidy;
  for(x=0;x<=(MESH_X-1);x++) {
    for (y=0; y<=(MESH_Y-1); y++) {
      gidy = x*MESH_Y + y;
      fscanf(fp, "%le %le %le\n",&x1,&y1,&gridinfo[gidy].compi[k]);
    }
  }
}
