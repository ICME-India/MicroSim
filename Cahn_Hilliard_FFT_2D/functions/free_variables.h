void free_variables(){

  long index, gidy, i;  
  long layer;
  
  
  //Free3M(Diffusivity, NUMPHASES, NUMCOMPONENTS-1);
  Free3M(ceq,         NUMPHASES, NUMPHASES);
  Free3M(cfill,       NUMPHASES, NUMPHASES);
  Free3M(AtomicMobility, NUMPHASES, NUMCOMPONENTS-1);

  fftw_free(phi);
  fftw_free(com);
  fftw_free(dfdphi);
  fftw_free(dfdc);
  fftw_destroy_plan(planF);
  fftw_destroy_plan(planB);

  index_count = layer_size*rows_x;
  for (index=0; index < index_count; index++) {
    free_memory_fields(&gridinfo[index]);
  }

  free(gridinfo);

  free(global_max_min.phi_max);
  free(global_max_min.phi_min);
  free(global_max_min.com_max);
  free(global_max_min.com_min);
  free(global_max_min.rel_change_phi);
  free(global_max_min.rel_change_com);

}
