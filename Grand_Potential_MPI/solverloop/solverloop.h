#ifndef SOLVERLOOP_H_
#define SOLVERLOOP_H_

void solverloop_phasefield_tdb(long *start, long *end);
void solverloop_concentration_tdb(long *start, long *end);
void collision_lbm(long *start, long *end);
void streaming_row(long *start, long *end);
void streaming_col(long *start, long *end);
void field_variables_lbm(long *start, long *end);

MPI_Request request_e1;
MPI_Request request_e2;
MPI_Request request_e3;
MPI_Request request_e4;
MPI_Request request_e5;
MPI_Request request_e6;
MPI_Request request_e7;
MPI_Request request_e8;
MPI_Request request_e9;
// void solverloop(long *start, long *end) {
//   long x;
//   
//   calculate_gradients(0, gradient);
//   swaplayers();
//   calculate_gradients(1, gradient);
//   swaplayers();
//   
//   for(x=start[X]-2; x<=end[X]+2; x++) {
//     calculate_gradients(x+1, gradient);    
//     calculate_divergence_phasefield(x, gradient);
//       
//   //Updating concentrations for the layer x-1
//     if (x > 1) {
//       calculate_fluxes_concentration(x-1, gradient);
//       calculate_divergence_concentration(x-1, gradient);
//     }
//     swaplayers();
//   }
// }
// void solverloop_phasefield_F_01(long *start, long *end) {
//   long x;
//   
//   calculate_gradients_phasefield(0, gradient, 1);
//   swaplayers();
//   calculate_gradients_phasefield(1, gradient, 1);
//   swaplayers();
//   
//   for(x=start[X]-2; x<=end[X]+2; x++) {
//     calculate_gradients_phasefield(x+1, gradient, 1);
//     calculate_divergence_phasefield(x, gradient);
//     swaplayers();
//   }
// }
// void solverloop_concentration_F_01(long *start, long *end) {
//   long b,k;
//   long x;
//   long INTERFACE_POS;
//   workers_max_min.INTERFACE_POS_MAX = 0;
//   
//   for (b=0; b < NUMPHASES; b++) {
//     workers_max_min.rel_change_phi[b] = 0.0;
//   }
//   for (k=0; k < NUMCOMPONENTS-1; k++) {
//     workers_max_min.rel_change_mu[k] = 0.0;
//   }
//   
//   calculate_gradients_concentration(0, gradient);
//   swaplayers();
//   calculate_gradients_concentration(1, gradient);
//   swaplayers();
//   
//   for(x=start[X]-2; x<=end[X]+2; x++) {
//     calculate_gradients_concentration(   x+1, gradient);
//     calculate_gradients_phasefield(      x+1, gradient);      
//   //Updating concentrations for the layer x-1
//     if (x > 1) {
//       calculate_fluxes_concentration(    x-1, gradient);
//       calculate_divergence_concentration(x-1, gradient);
// /***********************************************************************************/
//       if(SHIFT) {
//         //Check condition for the shift only for the lowest level
//         INTERFACE_POS = check_SHIFT(x-1);
//         if (INTERFACE_POS > workers_max_min.INTERFACE_POS_MAX) {
//           workers_max_min.INTERFACE_POS_MAX = INTERFACE_POS;
//         }
//       }
// /*********************************************************************************/
//     }
//     swaplayers();
//   }
// }

void collision_lbm(long *start, long *end) {
  long x;
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_distribution_function(x);
  }
}
void streaming_row_forward(long *start, long *end) {
  long x;
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_streaming_row_forward(x);
  }
}
void streaming_row_back(long *start, long *end) {
  long x;
  for(x=end[X]+2; x >= start[X]-2; x--) {
    calculate_streaming_row_back(x);
  }
}
void streaming_col_forward(long *start, long *end) {
  long x;
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_streaming_col_forward(x);
  }
}
void streaming_col_back(long *start, long *end) {
  long x;
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_streaming_col_back(x);
  }
}

void streaming_Z_forward(long *start, long *end) {
  long x;
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_streaming_Z_forward_3D(x);
  }
}

void streaming_Z_back(long *start, long *end) {
  long x;
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_streaming_Z_back_3D(x);
  }
}

void field_variables_lbm(long *start, long *end) {
  long x;
  for(x=start[X]-2; x<=end[X]+2; x++) {
    compute_field_variables_LBM(x);
  }
}


void iterative_stress_solver(long *start, long *end) {
  workers_avg_stiffness.C11 = 0.0;
  workers_avg_stiffness.C12 = 0.0;
  workers_avg_stiffness.C44 = 0.0;
  workers_avg_stress.xx   = 0.0;
  workers_avg_stress.yy   = 0.0;
  workers_avg_stress.xy   = 0.0;
  
  calculate_gradients_stress(0, gradient);
  swaplayers();
  calculate_gradients_stress(1, gradient);
  swaplayers();
  long x;
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_gradients_stress (x+1, gradient);
    calculate_divergence_stress(x,   gradient);
    swaplayers();
  }
}
void calculate_avg_strain() {
  MPI_Iallreduce(&workers_avg_stiffness.C11,  &avg_stiffness.C11,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request_e1);
  MPI_Iallreduce(&workers_avg_stiffness.C12,  &avg_stiffness.C12,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request_e2);
  MPI_Iallreduce(&workers_avg_stiffness.C44,  &avg_stiffness.C44,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request_e3);
  
  
  MPI_Iallreduce(&workers_avg_stress.xx,  &avg_stress.xx,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request_e4);
  MPI_Iallreduce(&workers_avg_stress.yy,  &avg_stress.yy,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request_e5);
  MPI_Iallreduce(&workers_avg_stress.xy,  &avg_stress.xy,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request_e6);
  
  MPI_Wait(&request_e1, MPI_STATUS_IGNORE);
  MPI_Wait(&request_e2, MPI_STATUS_IGNORE);
  MPI_Wait(&request_e3, MPI_STATUS_IGNORE);
  MPI_Wait(&request_e4, MPI_STATUS_IGNORE);
  MPI_Wait(&request_e5, MPI_STATUS_IGNORE);
  MPI_Wait(&request_e6, MPI_STATUS_IGNORE);
  
  if(DIMENSION == 3) {
    MPI_Iallreduce(&workers_avg_stress.zz,  &avg_stress.zz,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request_e7);
    MPI_Iallreduce(&workers_avg_stress.xz,  &avg_stress.xz,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request_e8);
    MPI_Iallreduce(&workers_avg_stress.yz,  &avg_stress.yz,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request_e9);
    
    MPI_Wait(&request_e7, MPI_STATUS_IGNORE);
    MPI_Wait(&request_e8, MPI_STATUS_IGNORE);
    MPI_Wait(&request_e9, MPI_STATUS_IGNORE);
  }
  
  
  if(DIMENSION == 2) {
    avg_stress.xx /= (MESH_X*MESH_Y);
    avg_stress.yy /= (MESH_X*MESH_Y);
    avg_stress.xy /= (MESH_X*MESH_Y);
    
    avg_stiffness.C11 /= (MESH_X*MESH_Y);
    avg_stiffness.C12 /= (MESH_X*MESH_Y);
    avg_stiffness.C44 /= (MESH_X*MESH_Y);
    
    if ((boundary[0][3].type == 2) || (boundary[1][3].type == 2)) {
      avg_stress.xx = boundary[0][3].value[0]/stiffness_phase[NUMPHASES-1].C44 - avg_stress.xx;
    }
    if ((boundary[2][3].type == 2) || (boundary[3][3].type == 2)) {
      avg_stress.yy = boundary[2][3].value[0]/stiffness_phase[NUMPHASES-1].C44 - avg_stress.yy;
    }
    avg_strain.xx = (avg_stiffness.C11*avg_stress.xx  - avg_stiffness.C12*avg_stress.yy)/(avg_stiffness.C11*avg_stiffness.C11 - avg_stiffness.C12*avg_stiffness.C12);
    avg_strain.yy = (-avg_stiffness.C12*avg_stress.xx + avg_stiffness.C11*avg_stress.yy)/(avg_stiffness.C11*avg_stiffness.C11 - avg_stiffness.C12*avg_stiffness.C12);
    
//     printf("avg_strain.yy=%le,avg_strain.xx=%le, avg_stress.yy=%le, avg_stress.xx=%le\n",avg_strain.yy,avg_strain.xx,avg_stress.yy, avg_stress.xx);
    
  } else {
    avg_stress.xx /= (MESH_X*MESH_Y*MESH_Z);
    avg_stress.yy /= (MESH_X*MESH_Y*MESH_Z);
    avg_stress.zz /= (MESH_X*MESH_Y*MESH_Z);
    avg_stress.xy /= (MESH_X*MESH_Y*MESH_Z);
    avg_stress.xz /= (MESH_X*MESH_Y*MESH_Z);
    avg_stress.yz /= (MESH_X*MESH_Y*MESH_Z);
    
    avg_stiffness.C11 /= (MESH_X*MESH_Y*MESH_Z);
    avg_stiffness.C12 /= (MESH_X*MESH_Y*MESH_Z);
    avg_stiffness.C44 /= (MESH_X*MESH_Y*MESH_Z);

    if ((boundary[0][3].type == 2) || (boundary[1][3].type == 2)) {
      avg_stress.xx = boundary[0][3].value[0]/stiffness_phase[NUMPHASES-1].C44 - avg_stress.xx;
    }
    if ((boundary[2][3].type == 2) || (boundary[3][3].type == 2)) {
      avg_stress.yy = boundary[2][3].value[0]/stiffness_phase[NUMPHASES-1].C44 - avg_stress.yy;
    }
    if ((boundary[4][3].type == 2) || (boundary[5][3].type == 2)) {
      avg_stress.zz = boundary[4][3].value[0]/stiffness_phase[NUMPHASES-1].C44 - avg_stress.zz;
    }
    /// newly added
    double det    = pow(avg_stiffness.C11,2) + avg_stiffness.C11 * avg_stiffness.C12 - 2 * pow(avg_stiffness.C12,2);

    avg_strain.xx =(+ avg_stress.xx * (avg_stiffness.C11 + avg_stiffness.C12) 
                    - avg_stress.yy * avg_stiffness.C12
                    - avg_stress.zz * avg_stiffness.C12) / det;

    avg_strain.yy =(- avg_stress.xx * avg_stiffness.C12 
                    + avg_stress.yy * (avg_stiffness.C11 + avg_stiffness.C12)
                    - avg_stress.zz * avg_stiffness.C12) / det;

    avg_strain.zz =(- avg_stress.xx * avg_stiffness.C12 
                    - avg_stress.yy * avg_stiffness.C12
                    + avg_stress.zz * (avg_stiffness.C11 + avg_stiffness.C12)) / det;
    
    avg_strain.xy = 0.5 * avg_stress.xy / avg_stiffness.C44;
    
    avg_strain.xz = 0.5 * avg_stress.xz / avg_stiffness.C44;

    avg_strain.yz = 0.5 * avg_stress.yz / avg_stiffness.C44;
  }
}

void solverloop_phasefield_tdb(long *start, long *end) {
  long x, b;
  long INTERFACE_POS;
  workers_max_min.INTERFACE_POS_MAX = 0;
  
  for (b=0; b < NUMPHASES; b++) {
    workers_max_min.rel_change_phi[b] = 0.0;
    workers_max_min.rel_change_phi_max[b] = 0.0;
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    workers_max_min.rel_change_mu[k] = 0.0;
    workers_max_min.rel_change_mu_max[k] = 0.0;
  }
  
  
//   printf("rank=%d, rank_x=%d, rank_y=%d, rank_z=%d\n", taskid, workers_mpi.rank_x, workers_mpi.rank_y, workers_mpi.rank_z);
  if ((FUNCTION_F !=5) && (!GRAIN_GROWTH)) {
    calculate_diffusion_potential( 0, gradient);
  }
  calculate_gradients_phasefield(0, gradient, 0);
  swaplayers();
  
  if ((FUNCTION_F !=5) && (!GRAIN_GROWTH)) {
    calculate_diffusion_potential( 1, gradient);
  }
  calculate_gradients_phasefield(1, gradient, 0);
  swaplayers();
  
  for(x=start[X]-2; x<=end[X]+2; x++) {
    if ((FUNCTION_F != 5) && (!GRAIN_GROWTH)) {
      calculate_diffusion_potential( x+1, gradient);
    }
    calculate_gradients_phasefield(x+1, gradient, 0);
    calculate_divergence_phasefield(x, gradient);
    if (FUNCTION_F == 5) {
      if(SHIFT) {
        //Check condition for the shift only for the lowest level
        INTERFACE_POS = check_SHIFT(x-1);
        if (INTERFACE_POS > workers_max_min.INTERFACE_POS_MAX) {
          workers_max_min.INTERFACE_POS_MAX = INTERFACE_POS;
        }
      }
    }
    swaplayers();
  }
}
void solverloop_concentration_tdb(long *start, long *end) {
  long b,k;
  long x;
  long INTERFACE_POS;
  workers_max_min.INTERFACE_POS_MAX = 0;
  
  for (b=0; b < NUMPHASES; b++) {
    workers_max_min.rel_change_phi[b] = 0.0;
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    workers_max_min.rel_change_mu[k] = 0.0;
  }

  
//   for (b=0; b < NUMPHASES; b++) {
//     workers_max_min.rel_change_phi[b] = 0.0;
//   }
//   for (k=0; k < NUMCOMPONENTS-1; k++) {
//     workers_max_min.rel_change_mu[k] = 0.0;
//   }
  
  calculate_gradients_phasefield(   0, gradient, 1); 
  swaplayers();
  calculate_gradients_phasefield(   1, gradient, 1); 
  calculate_gradients_concentration(0, gradient);
//   calculate_gradients_concentration(1, gradient);
  swaplayers();
  
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_gradients_phasefield(      x+1, gradient, 1); 
    calculate_gradients_concentration(   x,   gradient);
  //Updating concentrations for the layer x-1
    if (x > 1) {
      calculate_fluxes_concentration(    x-1, gradient);
      calculate_divergence_concentration(x-1, gradient);
/***********************************************************************************/
      if(SHIFT) {
        //Check condition for the shift only for the lowest level
        INTERFACE_POS = check_SHIFT(x-1);
        if (INTERFACE_POS > workers_max_min.INTERFACE_POS_MAX) {
          workers_max_min.INTERFACE_POS_MAX = INTERFACE_POS;
        }
      }
/*********************************************************************************/
    }
    swaplayers();
  }
}
void  smooth(long *start, long *end) {
  long x;
//   calculate_diffusion_potential( 0, gradient);
  calculate_gradients_phasefield(   0, gradient, 1); 
  swaplayers();
//   calculate_diffusion_potential( 1, gradient);
  calculate_gradients_phasefield(   1, gradient, 1);
  if ((FUNCTION_F!=5) && (!GRAIN_GROWTH)) {
    calculate_gradients_concentration(0, gradient);
  }
//   calculate_gradients_concentration(1, gradient);
  swaplayers();
  
  for(x=start[X]-2; x<=end[X]+2; x++) {
//     calculate_gradients(x+1, gradient);
//     calculate_divergence_phasefield_smooth(x, gradient);
//     calculate_diffusion_potential( x+1, gradient);
    calculate_gradients_phasefield(        x+1, gradient, 1);
    if ((FUNCTION_F !=5) && (!GRAIN_GROWTH)) {
      calculate_gradients_concentration(     x,   gradient);
    }
    calculate_divergence_phasefield_smooth(x,   gradient);  
  //Updating concentrations for the layer x-1
    if (x > 1) {
      if ((FUNCTION_F !=5) && (!GRAIN_GROWTH)) {
        calculate_divergence_concentration_smooth(x-1, gradient);
      }
    }
    swaplayers();
  }
}
// void  smooth_concentration(long *start, long *end) {
//   long x;
//   
//   calculate_gradients(0, gradient);
//   swaplayers();
//   calculate_gradients(1, gradient);
//   swaplayers();
//   
//   for(x=start[X]-2; x<=end[X]+2; x++) {
//     calculate_gradients(x+1, gradient);      
//   //Updating concentrations for the layer x-1
//     if (x > 1) {
//       calculate_divergence_concentration_smooth_concentration(x-1, gradient);
//     }
//     swaplayers();
//   }
// }
#endif
