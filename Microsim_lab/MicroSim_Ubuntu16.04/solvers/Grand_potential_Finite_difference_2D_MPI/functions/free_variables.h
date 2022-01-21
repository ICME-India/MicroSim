void free_variables(){

  long index, gidy, i;  
  long index_count;
  long layer;
  
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
  FreeM(Gamma, NUMPHASES);
  Free3M(Gamma_abc, NUMPHASES, NUMPHASES);
  
  if(FUNCTION_ANISOTROPY !=0) {
    if(FOLD==4) {
      FreeM(dab, NUMPHASES);
    }
  }
  
  free(C);

  for(layer=0; layer < 4; layer++) {
    for (gidy=0; gidy < workers_mpi.layer_size; gidy++) {
      free_memory_gradlayer(&gradient1[layer][gidy]);
    }
    free(gradient1[layer]);
  }
  
  free_memory_fields(gridinfo_instance);
  free(gridinfo_instance);
  
  index_count = workers_mpi.layer_size*workers_mpi.rows_x;
  
  for (index=0; index < index_count; index++) {
    if ((&gridinfo_w[index]) !=NULL) {
      free_memory_fields(&gridinfo_w[index]);
    }
  }

  free(gridinfo_w);
  free(eigen_strain);
  free(Stiffness_c);
  free(Stiffness_t);
  
  
  for (i=0; i<6; i++) {
    free(boundary[i]);
  }
  
  if (taskid == MASTER) {
    free(global_max_min.phi_max);
    free(global_max_min.phi_min);
    free(global_max_min.mu_max);
    free(global_max_min.mu_min);
    free(global_max_min.rel_change_phi);
    free(global_max_min.rel_change_mu);
  }
  
  free(workers_max_min.phi_max);
  free(workers_max_min.phi_min);
  free(workers_max_min.mu_max);
  free(workers_max_min.mu_min);
  free(workers_max_min.rel_change_phi);
  free(workers_max_min.rel_change_mu);
  
  

  free(buffer_boundary_x);
  free(buffer_boundary_y);

  
  Free3M(cmu,NUMCOMPONENTS-1,NUMCOMPONENTS-1);
  Free3M(muc,NUMCOMPONENTS-1,NUMCOMPONENTS-1);
//   
  FreeM(dcdmu,    NUMCOMPONENTS-1);
  FreeM(inv_dcdmu,NUMCOMPONENTS-1);
  Free4M(Rotation_matrix,NUMPHASES, NUMPHASES, DIMENSION);
  Free4M(Inv_Rotation_matrix,NUMPHASES, NUMPHASES, DIMENSION);
  free(Rotated_qab);
  free(deltamu);
  free(sum);
  free(start);
  free(end);
  free(averow);
  free(rows);
  free(offset);
  free(extra);
  free(divphi);
  free(lambda_phi);
  free(divflux);
  free(divjat);
  free(deltac);
}