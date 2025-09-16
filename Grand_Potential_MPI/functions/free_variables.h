
/// @brief Frees allocated memory
void free_variables()
{
  long index, gidy, i;  
  long index_count;
  long layer;

  if (USE_NEW_INPFILE)
  {
    ms_GlobalFree_PhaseFieldMatrices();
  }
  else 
  {
    free(start);//
    free(end);//
    free(averow);//
    free(rows);//
    free(offset);//
    free(extra);//

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
    FreeM(Gamma,        NUMPHASES);

    if (NUMPHASES > 2)
    {
      Free3M(Gamma_abc,   NUMPHASES, NUMPHASES);      
    }
    
    Free3M(c_guess,     NUMPHASES, NUMPHASES);    

    Free4M(Rotation_matrix,NUMPHASES, NUMPHASES, DIMENSION);
    Free4M(Inv_Rotation_matrix,NUMPHASES, NUMPHASES, DIMENSION);
    
    free(Rotated_qab);

    if(dab) {
      FreeM(dab, NUMPHASES);
    }

    if(USING_FAB) {
      FreeM(fab, NUMPHASES);
    }
    
    free(C);

    Free3M(cmu,NUMPHASES,NUMCOMPONENTS-1);
    Free3M(muc,NUMPHASES,NUMCOMPONENTS-1);

    for (i = 0; i < NUMPHASES; ++i) {
      free(Phases[i]);
    }
    free(Phases);
    Phases = NULL;
    
    for (i = 0; i < NUMCOMPONENTS; ++i) {
      free(Components[i]);
    }
    free(Components);
    Components = NULL;
    
    for (i = 0; i < NUM_THERMO_PHASES; ++i) {
      free(Phases_tdb[i]);
    }
    free(Phases_tdb);
    Phases_tdb = NULL;
    
    for (i = 0; i < NUMPHASES; ++i) {
      free(phase_map[i]);
    }
    free(phase_map);
    phase_map = NULL;
    
    free(thermo_phase);
  }

  if (FUNCTION_F==4) {
    for(i=0;i<NUM_THERMO_PHASES;i++) {
      for(j=0;j<NUMCOMPONENTS-1;j++){
        for(k=0;k<NUMCOMPONENTS-1;k++){
          gsl_spline_free (spline_ThF[i][j][k]);
          gsl_interp_accel_free (acc_ThF[i][j][k]);
        }
        free(spline_ThF[i][j]);
        free(acc_ThF[i][j]);
      }
//       free(Mat[i]);
      free(spline_ThF[i]);
      free(acc_ThF[i]);
    }
    free(spline_ThF);
    free(acc_ThF);
//     Mat=NULL;
    
    for(i=0;i<NUM_THERMO_PHASES-1;i++) {
      for(j=0;j<NUMCOMPONENTS-1;j++){
        for(k=0;k<2;k++){
          gsl_spline_free (spline_ES[i][j][k]);
          gsl_interp_accel_free (acc_ES[i][j][k]);
        }
        free(spline_ES[i][j]);
        free(acc_ES[i][j]);
      }
//       free(Mat[i]);
      free(spline_ES[i]);
      free(acc_ES[i]);
    }
    free(spline_ES);
    free(acc_ES);
  }

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

  // exit(0);

  free(gridinfo_w);
  free(eigen_strain_phase);
  free(stiffness_phase);
  free(stiffness_phase_n);
  free(stiffness_t_phase);
  
  if (ELASTICITY) {
    free(iter_gridinfo_w);
  }
  
  if (LBM) {
    for (index=0; index < index_count; index++) {
      if ((&lbm_gridinfo_w[index]) !=NULL) {
        free_memory_lbm_fields(&lbm_gridinfo_w[index]);
      }
//       free(temp_fdis[index]);
    }
//     free(lbm_gridinfo_w);
    free(rho_LBM);
    free(beta_c);
//     free(temp_fdis);
    Free3M(gradc,workers_mpi.layer_size,NUMCOMPONENTS-1);
    free(lbm_gridinfo_w);
  }
  
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
  free(buffer_boundary_z);
  
  free(buffer_boundary_x_stress);
  free(buffer_boundary_y_stress);
  free(buffer_boundary_z_stress);
  
  Free3M(dcdmu_phase, NUMPHASES, NUMCOMPONENTS-1);
 
//   
  FreeM(dcdmu,    NUMCOMPONENTS-1);
  FreeM(inv_dcdmu,NUMCOMPONENTS-1);
  FreeM(Ddcdmu,   NUMCOMPONENTS-1);
  
  free(deltamu);
  free(sum);
  
  free(divphi);
  free(lambda_phi);
  free(divflux);
  free(divjat);
  free(deltac);
  free(c_old);
  free(c_new);
  free(c_tdt);
  
}
