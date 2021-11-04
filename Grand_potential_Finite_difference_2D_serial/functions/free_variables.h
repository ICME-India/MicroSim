void free_variables(){

  long index, gidy, i;  
  long index_count;
  long layer;
  
  index_count = layer_size*rows_x;
  
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

//   
//   free_memory_gradlayer(grad);
//   free(grad);
  free_memory_fields(gridinfo_instance);
  for (index=0; index < index_count; index++) {
    free_memory_fields(&gridinfo[index]);
  }
  
  for(layer=0; layer < 4; layer++) {
    for (gidy=0; gidy < layer_size; gidy++) {
      free_memory_gradlayer(&gradient1[layer][gidy]);
    }
    free(gradient1[layer]);
  }
  free(gridinfo);
  free(eigen_strain);
  free(Stiffness_c);
  free(Stiffness_t);
//   
  for (i=0; i<6; i++) {
    free(boundary[i]);
  }
  
  free(global_max_min.phi_max);
  free(global_max_min.phi_min);
  free(global_max_min.mu_max);
  free(global_max_min.mu_min);
  free(global_max_min.rel_change_phi);
  free(global_max_min.rel_change_mu);
  
  Free3M(cmu,NUMCOMPONENTS-1,NUMCOMPONENTS-1);
  Free3M(muc,NUMCOMPONENTS-1,NUMCOMPONENTS-1);
  
  FreeM(dcdmu,    NUMCOMPONENTS-1);
  FreeM(inv_dcdmu,NUMCOMPONENTS-1);
  Free4M(Rotation_matrix,NUMPHASES, NUMPHASES, DIMENSION);
  Free4M(Inv_Rotation_matrix,NUMPHASES, NUMPHASES, DIMENSION);
  free(Rotated_qab);
  free(deltamu);
  free(deltac);
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
}