void GEwriter_Execution(struct variables* gridinfo,int t) {
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
