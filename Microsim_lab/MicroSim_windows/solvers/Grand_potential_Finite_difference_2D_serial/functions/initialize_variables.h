#ifndef INITIALIAZE_VARIABLES_H_
#define INITIALIAZE_VARIABLES_H_

void initialize_variables() {
  long a, i, j;
  for (a=0;a<NUMPHASES;a++) {
    for (i=0;i<NUMCOMPONENTS-1;i++) {
      for (j=0;j<NUMCOMPONENTS-1;j++) {
	if (i==j) {
	  muc[a][i][j]=2.0*A[a][i][j];
	} else {
	  muc[a][i][j]=A[a][i][j];
	}
      }
    }
    matinvnew(muc[a], cmu[a], NUMCOMPONENTS-1);
  }
  
// #ifndef ISOTHERMAL
  for (a=0;a<NUMPHASES;a++) {
    for (i=0; i < NUMCOMPONENTS-1; i++) {      
      dcbdT_phase[a][i] = 0.0;
      for (j=0; j < NUMCOMPONENTS-1; j++) {
        dcbdT_phase[a][i] += cmu[a][i][j]*(-dBbdT[a][j]);
//           dcbdT_phase[a][i] += dc_dmu(a,i,j)*(-dBbdT[a][j]);
      }
    }
  }
// #endif
  dcdmu      = MallocM((NUMCOMPONENTS-1),(NUMCOMPONENTS-1));
  inv_dcdmu  = MallocM(NUMCOMPONENTS-1,NUMCOMPONENTS-1);
  deltamu    = MallocV(NUMCOMPONENTS-1);
  deltac     = MallocV(NUMCOMPONENTS-1);
  sum        = MallocV(NUMCOMPONENTS-1);
  divphi     = MallocV(NUMPHASES);
  lambda_phi = MallocV(NUMPHASES);
  divflux    = MallocV(NUMCOMPONENTS-1);
  c_old      = MallocV(NUMCOMPONENTS-1);
  c_new      = MallocV(NUMCOMPONENTS-1);
  divjat     = MallocV(NUMCOMPONENTS-1);
  deltac     = MallocV(NUMCOMPONENTS-1);
  
  //Need to initialize according to the dimension. 2D rotation matrices and 3D rotation matrices
  //Need to take in values of the orientation of the solid phases with a common reference frame
  
//   for (a=0; a < NUMPHASES; a++) {
// //           Rth_phase                                 =  Rtheta*(double)(rand())/(double)RAND_MAX;
//     Rth_phase                                 =  a*fabs(Rtheta) + Rtheta;
//     Rotation_matrix[a][NUMPHASES-1][0][0]     =  cos(Rth_phase*M_PI/180.0);
//     Rotation_matrix[a][NUMPHASES-1][0][1]     = -sin(Rth_phase*M_PI/180.0);
//     Rotation_matrix[a][NUMPHASES-1][1][0]     =  sin(Rth_phase*M_PI/180.0);
//     Rotation_matrix[a][NUMPHASES-1][1][1]     =  cos(Rth_phase*M_PI/180.0);
//     
//     Rotation_matrix[NUMPHASES-1][a][0][0]     =  cos(Rth_phase*M_PI/180.0);
//     Rotation_matrix[NUMPHASES-1][a][0][1]     = -sin(Rth_phase*M_PI/180.0);
//     Rotation_matrix[NUMPHASES-1][a][1][0]     =  sin(Rth_phase*M_PI/180.0);
//     Rotation_matrix[NUMPHASES-1][a][1][1]     =  cos(Rth_phase*M_PI/180.0);
//     
//     
//     Inv_Rotation_matrix[a][NUMPHASES-1][0][0] =   cos(Rth_phase*M_PI/180.0);
//     Inv_Rotation_matrix[a][NUMPHASES-1][0][1] =   sin(Rth_phase*M_PI/180.0);
//     Inv_Rotation_matrix[a][NUMPHASES-1][1][0] =  -sin(Rth_phase*M_PI/180.0);
//     Inv_Rotation_matrix[a][NUMPHASES-1][1][1] =   cos(Rth_phase*M_PI/180.0);
//     
//     Inv_Rotation_matrix[NUMPHASES-1][a][0][0] =   cos(Rth_phase*M_PI/180.0);
//     Inv_Rotation_matrix[NUMPHASES-1][a][0][1] =   sin(Rth_phase*M_PI/180.0);
//     Inv_Rotation_matrix[NUMPHASES-1][a][1][0] =  -sin(Rth_phase*M_PI/180.0);
//     Inv_Rotation_matrix[NUMPHASES-1][a][1][1] =   cos(Rth_phase*M_PI/180.0);
//   }
//   for (a=0; a < NUMPHASES-1; a++) {
//     for (b=0; b < NUMPHASES-1; b++) {
//       Rotation_matrix[a][b][0][0]     =  0.0;
//       Rotation_matrix[a][b][0][1]     =  0.0;
//       Rotation_matrix[a][b][1][0]     =  0.0;
//       Rotation_matrix[a][b][1][1]     =  0.0;
//       
//       Inv_Rotation_matrix[a][b][0][0] =  0.0;
//       Inv_Rotation_matrix[a][b][0][1] =  0.0;
//       Inv_Rotation_matrix[a][b][1][0] =  0.0;
//       Inv_Rotation_matrix[a][b][1][1] =  0.0;
//     }
//   }
  
//   workers_mpi.lastx=0;
//   workers_mpi.firstx=0;
//   workers_mpi.lasty=0;
//   workers_mpi.firsty=0;
}

#endif
